import streamlit as st
import pandas as pd
from itertools import permutations
import numpy as np
import io
import base64
from io import BytesIO
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from concurrent.futures import ProcessPoolExecutor
import os

# Função para download de DataFrame como xlsx
def download_xlsx(df, filename):
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name='Sheet1', index=False)
    output.seek(0)
    st.download_button(
        label="Download Data as XLSX",
        data=output,
        file_name=filename,
        mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
    )

# Função para calcular os índices DRIS e IBN
def calculate_dris_ibn(df, nutrients, response_variable):
    pairs = list(permutations(nutrients, 2))
    all_ratios = {}
    for pair in pairs:
        ratio_name = f'{pair[0]}/{pair[1]}'
        all_ratios[ratio_name] = df[pair[0]] / df[pair[1]]
    df_ratios = pd.concat([df, pd.DataFrame(all_ratios)], axis=1)
    threshold = df[response_variable].mean() + 0.5 * df[response_variable].std()
    df_high_yield = df[df[response_variable] >= threshold]
    df_low_yield = df[df[response_variable] < threshold]
    df_high_yield_ratios = df_high_yield.copy()
    df_low_yield_ratios = df_low_yield.copy()
    for pair in pairs:
        df_high_yield_ratios[f"{pair[0]}/{pair[1]}"] = df_high_yield_ratios[pair[0]] / df_high_yield_ratios[pair[1]]
        df_low_yield_ratios[f"{pair[0]}/{pair[1]}"] = df_low_yield_ratios[pair[0]] / df_low_yield_ratios[pair[1]]
    var_high_yield_ratios = df_high_yield_ratios.iloc[:, -len(pairs):].var()
    var_low_yield_ratios = df_low_yield_ratios.iloc[:, -len(pairs):].var()
    variance_ratio_ratios = var_low_yield_ratios / var_high_yield_ratios
    max_variance_ratio_ratios = {}
    all_variance_ratio_ratios = {nutrient: {} for nutrient in nutrients}
    for nutrient in nutrients:
        variance_ratios_current_nutrient = {}
        for other_nutrient in nutrients:
            if nutrient != other_nutrient:
                variance_ratio_1 = variance_ratio_ratios.get(f'{nutrient}/{other_nutrient}', 0)
                variance_ratio_2 = variance_ratio_ratios.get(f'{other_nutrient}/{nutrient}', 0)
                if variance_ratio_1 > variance_ratio_2:
                    variance_ratios_current_nutrient[f'{nutrient}/{other_nutrient}'] = variance_ratio_1
                else:
                    variance_ratios_current_nutrient[f'{other_nutrient}/{nutrient}'] = variance_ratio_2
        max_variance_ratio_ratios[nutrient] = max(variance_ratios_current_nutrient, key=variance_ratios_current_nutrient.get)
        all_variance_ratio_ratios[nutrient] = variance_ratios_current_nutrient
    dris_indices = {}
    for nutrient in nutrients:
        f_high_total = 0
        f_low_total = 0
        for ratio, variance_ratio in all_variance_ratio_ratios[nutrient].items():
            numer, denom = ratio.split('/')
            if numer == nutrient:
                Ea_num = (df_low_yield[numer] / df_low_yield[denom]).mean()
                En_num = (df_high_yield[numer] / df_high_yield[denom]).mean()
                s_num = (df_high_yield[numer] / df_high_yield[denom]).std()
                f_high = f_E_Ei(Ea_num, En_num, s_num)
                f_high_total += f_high
            else:
                Ea_den = (df_low_yield[numer] / df_low_yield[nutrient]).mean()
                En_den = (df_high_yield[numer] / df_high_yield[nutrient]).mean()
                s_den = (df_high_yield[numer] / df_high_yield[nutrient]).std()
                f_low = f_E_Ei(Ea_den, En_den, s_den)
                f_low_total += f_low
        dris_index = (f_high_total - f_low_total) / (len(nutrients) - 1)
        dris_indices[nutrient] = dris_indices.get(nutrient, 0) + dris_index
    IBN = sum(abs(value) for value in dris_indices.values()) / len(nutrients)
    return dris_indices, IBN, df_ratios, df_high_yield, all_variance_ratio_ratios

def f_E_Ei(E_Ea, E_En, s, K=10):
    return (E_Ea - E_En) * (K / s)

def bootstrap_iteration(df, nutrients, response_variable):
    df_sample = df.sample(n=len(df), replace=True)
    dris_indices_sample, IBN_sample, _, _, _ = calculate_dris_ibn(df_sample, nutrients, response_variable)
    return dris_indices_sample, IBN_sample

def process_tab2(df, nutrients, response_variable, bootstrap_iterations):
    with ProcessPoolExecutor(max_workers=os.cpu_count() - 1) as executor:
        futures = [executor.submit(bootstrap_iteration, df, nutrients, response_variable) for _ in range(bootstrap_iterations)]
        results = [future.result() for future in futures]
    return results

def bootstrap_iteration_tab5(df_high_yield, dris_indices_high_yield_df, nutrient):
    df_sample = df_high_yield.sample(n=len(df_high_yield), replace=True)
    dris_sample = pd.DataFrame({nutrient: dris_indices_high_yield_df[nutrient].sample(n=len(dris_indices_high_yield_df), replace=True)})
    coeffs = np.polyfit(dris_sample[nutrient], df_sample[nutrient], 1)
    zero = coeffs[1]
    return zero

def process_tab5(df_high_yield, dris_indices_high_yield_df, nutrients, bootstrap_iterations):
    results = {}
    for nutrient in nutrients:
        with ProcessPoolExecutor(max_workers=os.cpu_count() - 1) as executor:
            futures = [executor.submit(bootstrap_iteration_tab5, df_high_yield, dris_indices_high_yield_df, nutrient) for _ in range(bootstrap_iterations)]
            zeros_samples = [future.result() for future in futures]
        results[nutrient] = zeros_samples
    df_zeros_samples = pd.DataFrame(results)
    return df_zeros_samples


def plot_kde(nutrient, samples, pdf=None):
    plt.figure(figsize=(8, 4))
    sns.kdeplot(samples[nutrient], color='green', fill=True)
    plt.xlabel(f'{nutrient}')
    plt.ylabel(translate("density_label"))
    plt.title(translate("bootkde_label") + f' {nutrient}')
    lower_bound = np.percentile(samples[nutrient], 2.5)
    upper_bound = np.percentile(samples[nutrient], 97.5)
    plt.axvline(lower_bound, color='red', linestyle='--', linewidth=1)
    plt.axvline(upper_bound, color='red', linestyle='--', linewidth=1)
    st.pyplot(plt.gcf())

    if pdf:
        pdf.savefig(plt.gcf())
    
    pdf_buffer = io.BytesIO()
    plt.savefig(pdf_buffer, format='pdf')
    st.download_button(label=translate("dowlkde_label") + " " + nutrient + " " + translate("as_PDF_label"),
                       data=pdf_buffer.getvalue(),
                       file_name=f"{nutrient}_bootstrap_plot.pdf",
                       mime="application/pdf")

    svg_buffer = io.BytesIO()
    plt.savefig(svg_buffer, format='svg')
    st.download_button(label=translate("dowlkde_label") + " " + f'{nutrient}' + " " + translate("as_SVG_label"),
                       data=svg_buffer.getvalue(),
                       file_name=f"{nutrient}_bootstrap_plot.svg",
                       mime="image/svg+xml")

    plt.close()
    
def plot_original_regression(nutrient):
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x=dris_indices_high_yield_df[nutrient], y=df_high_yield[nutrient], color='blue')
    sns.regplot(x=dris_indices_high_yield_df[nutrient], y=df_high_yield[nutrient], scatter=False, line_kws={'color': 'red'})
    plt.xlabel(f"DRIS-{nutrient}")
    plt.ylabel(f"{nutrient}")
    slope, intercept = np.polyfit(dris_indices_high_yield_df[nutrient], df_high_yield[nutrient], 1)
    r_squared = np.corrcoef(dris_indices_high_yield_df[nutrient], df_high_yield[nutrient])[0, 1] ** 2
    import matplotlib
    from matplotlib import mathtext
    matplotlib.rc('text', usetex=False)
    matplotlib.rc('font', family='serif')
    matplotlib.rc('mathtext', fontset='cm')
    math_expr = r'$y^{\hat{}}=' + rf'{intercept:.2f}+{slope:.2f}x$'
    plt.text(0.05, 0.95, math_expr, fontsize=12, transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=0.8))
    plt.text(0.05, 0.90, f"R² = {r_squared:.2f}", fontsize=12, transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=0.8))
    plt.tight_layout()
    st.pyplot(plt.gcf())
    regression_pdf.savefig(plt.gcf())
    regression_svg_buffer = io.BytesIO()
    plt.savefig(regression_svg_buffer, format='svg')
    st.download_button(label=translate("dowlreg_label") +" "+ f' {nutrient}'+" "+translate("as_SVG_label"), key=f"svg_dwltab51_{nutrient}_original",
                       data=regression_svg_buffer.getvalue(),
                       file_name=f"regression_plot_{nutrient}_original.svg",
                       mime="image/svg+xml")
    plt.close()
def add_bg_from_local(image_file):
    with open(image_file, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read()).decode()
    return f"data:image/png;base64,{encoded_string}"    
    

with open('data/fertilizer_8809685.png', 'rb') as file:
    image_bytes = file.read()
    base64_icon = base64.b64encode(image_bytes).decode('utf-8')





# Configuração da página
st.set_page_config(layout="wide", page_title="Análise DRIS", page_icon=f"data:image/png;base64,{base64_icon}")

# Carregar imagens de fundo
main_bg = add_bg_from_local('data/page2.jpeg')  # Substitua pelo caminho real
sidebar_bg = add_bg_from_local('data/sidebar2.jpg')  # Substitua pelo caminho real


# Adicionar CSS personalizado para as imagens de fundo
st.markdown(
    f"""
    <style>
    .main .block-container {{
        background-image: url("{main_bg}");
        background-size: cover;
        background-position: center;
        background-repeat: no-repeat;
    }}
    .main .block-container::before {{
        content: "";
        position: absolute;
        top: 0;
        left: 0;
        width: 100%;
        height: 100%;
        background-color: rgba(255, 255, 255, 0.5);  # Ajuste o último valor (0.5) para mudar a opacidade
        pointer-events: none;  # Permite que os elementos abaixo sejam clicáveis
    }}
    [data-testid="stSidebar"] {{
        background-image: url("{sidebar_bg}");
        background-size: cover;
        background-position: center;
        background-repeat: no-repeat;
    }}
    [data-testid="stSidebar"]::before {{
        content: "";
        position: absolute;
        top: 0;
        left: 0;
        width: 100%;
        height: 100%;
        background-color: rgba(255, 255, 255, 0.5);  # Ajuste o último valor (0.5) para mudar a opacidade
        pointer-events: none;  # Permite que os elementos abaixo sejam clicáveis
    }}
    .main-title {{
        font-size: 36px;
        color: #1E1E1E;
        font-weight: bold;
        text-align: center;
        margin-bottom: 20px;
    }}
    
    /* CSS para fonte preta, tamanho 20px e negrito */
    .stTextInput > div > div > input,
    .stSelectbox > div > div > select,
    .stMultiSelect > div > div > select,
    .stSlider > div > div > div > div,
    .stNumberInput > div > div > input,
    .stDateInput > div > div > input,
    .stTimeInput > div > div > input,
    .stTextArea > div > div > textarea {{
        color: black !important;
        font-size: 30px !important;
        font-weight: bold !important;
    }}

    /* Estilo para botões */
    .stButton > button {{
        color: black !important;
        font-size: 30px !important;
        font-weight: bold !important;
    }}

    /* Estilo para texto em geral */
    .stMarkdown, .stText {{
        color: black !important;
        font-size: 30px !important;
        font-weight: bold !important;
    }}

    /* Estilo para as abas */
    .stTabs [data-baseweb="tab"] {{
        color: #2C3E50 !important;
        font-size: 18px !important;
        
    }}
    .stTabs [data-baseweb="tab-list"] {{
        gap: 8px;
    }}
    .stTabs [data-baseweb="tab"] [data-testid="stMarkdownContainer"] p {{
        font-size: 18px !important;
        
    }}
    /* Estilo para o conteúdo dentro das abas */
    .stTabs [data-baseweb="tab-panel"] {{
        font-size: 18px !important;
        
        color: #2C3E50 !important;
    }}

        
  
    /* Estilo atualizado para o texto da barra lateral */
    [data-testid="stSidebar"] [data-testid="stMarkdown"] {{
        color: #330033 !important ;       
        font-size: 20px ;
        font-weight: bold !important;
    }}
    
    /* Estilos adicionais para garantir que o texto da barra lateral seja formatado corretamente */
    [data-testid="stSidebar"] [data-testid="stMarkdown"] p {{
        color: #330033 !important;
        font-size: 20px !important;
        font-weig    }}
    
    /* Estilo para links na barra lateral */
    [data-testid="stSidebar"] [data-testid="stMarkdown"] a {{
        color: #330033 !important;
        font-size: 20px !important;
        font-weight: bold !important;
    }}
    
    /* Estilo para cabeçalhos na barra lateral */
    [data-testid="stSidebar"] [data-testid="stMarkdown"] h1,
    [data-testid="stSidebar"] [data-testid="stMarkdown"] h2,
    [data-testid="stSidebar"] [data-testid="stMarkdown"] h3,
    [data-testid="stSidebar"] [data-testid="stMarkdown"] h4,
    [data-testid="stSidebar    [data-testid="stSidebar"] [data-testid="stMarkdown"] h6 {{
        color: #330033 !important;
        font-size: 20px !important;
        font-weight: bold !important;
    }}

    /* Estilo para listas na barra lateral */
    [data-testid="stSidebar"] [data-testid="s    [data-testid="stSidebar"] [data-testid="stMarkdown"] ol {{
        color: #330033 !important;
        font-size: 20px !important;
        font-weight: bold !important;
    }}
    
    </style>
    """,
    unsafe_allow_html=True
)

# Título principal
st.markdown('<h1 class="main-title">DRIS</h1>', unsafe_allow_html=True)

# Conteúdo da barra lateral
st.sidebar.markdown("**Versão:** 0.1")
st.sidebar.markdown("**Autor:** Prof. Walter E. Pereira, UFPB, CCA")

# Selecionar o idioma
language = st.sidebar.selectbox("Select Language", ["English", "Português", "Español"])

# Dicionário com traduções
translations = {
    "English": {
        "title": "DRIS Analysis",
        "upload_label": "Upload your data file (xlsx or csv)",
        "show_example_label": "Show Example File",
        "help_label": "Help",
        "response_variable_label": "Select the response variable (productivity):",
        "download_results_label": "Download Results",
        "dris_indices_label": "DRIS Indices",
        "ibn_label": "NBI",
        "means_and_sds_label": "Means and Standard Deviations for High Productivity",
        "dris_indices_high_yield_label": "DRIS Indices for High Productivity",
        "reference_values_label": "Reference Values",
        "bootstrap_iterations_label": "Number of Bootstrap iterations",
        "run_bootstrap_label": "Run Bootstrap",
        "example_file_not_found_error": "The example file '{example_file}' was not found.",
        "faqs_title": "Frequently Asked Questions (FAQs)",
        "faqs_q1": "What types of files can I upload?",
        "faqs_a1": "You can upload Excel (.xlsx) or CSV (.csv) files containing nutrient data and the response variable.",
        "faqs_q2": "How do I select the response variable?",
        "faqs_a2": "After uploading the file, choose the column that represents the response variable (productivity) in the 'Select the response variable (productivity):' dropdown menu.",
        "faqs_q3": "What are the DRIS indices?",
        "faqs_a3": "The DRIS indices measure the relationship between the concentration of different nutrients in the plant, indicating if the plant is well-nourished or if there is any imbalance.",
        "faqs_q4": "What is the NBI?",
        "faqs_a4": "The NBI (Nutritional Balance Index) is a measure that summarizes the degree of nutritional imbalance in the plant. A high NBI indicates a greater imbalance.",
        "faqs_q5": "What is bootstrap?",
        "faqs_a5": "Bootstrap analysis is a statistical technique that allows you to calculate confidence intervals for the DRIS indices and NBI.",
        "faqs_q6": "How do I interpret the reference values?",
        "faqs_a6": "The reference values indicate the ideal concentration of each nutrient in the plant to achieve high productivity.",
        "example_usage_title": "Detailed Example of Usage",
        "example_usage_step1": "1. Load the data file:",
        "example_usage_step1_download": "Download the example file [link for example file].",
        "example_usage_step1_upload": "Load the file into the application using the 'Upload your data file' button.",
        "example_usage_step2": "2. Select the response variable:",
        "example_usage_step2_select": "Select the 'Productivity' column as the response variable.",
        "example_usage_step3": "3. Run the analysis:",
        "example_usage_step3_run": "Click the 'Run Analysis' button to calculate the DRIS indices, NBI, and other relevant information.",
        "example_usage_step4": "4. View the results:",
        "example_usage_step4_view": "Explore the tabs 'DRIS Indices', 'Bootstrap DRIS/NBI', 'Means and Standard Deviations', 'DRIS Indices for High Productivity', and 'Reference Values' to see the analysis results.",
        "example_usage_step5": "5. Download the results:",
        "example_usage_step5_download": "Use the download buttons to save the results to Excel files.",
        "select_language_label": "Select Language",
        "density_label": "Density",
        "bootkde_label": "Bootstrap distribution of ",
        "dowlkde_label": "Download Bootstrap plot of ",
        "dowlreg_label": "Download regression plot of ",
        "as_PDF_label": "as PDF",
        "as_SVG_label": "as SVG",
        "dldrisridge_label": "Download DRIS Ridge Plot as PDF",
        "dldrisridgesvg_label": "Download DRIS Ridge Plot as SVG",
        "reference_label": "Download All Reference Value Plots as PDF",
        "allreg_label": "Download All Regression Plots as PDF",
        "dwlibnpdf_label": "Download NBI Plot as PDF",
        "dwlibnsvg_label": "Download NBI Plot as SVG",
        "dris_indices_label1": "summary of DRIS Indices",
        "dris_indices_label2": "full data of DRIS Indices",
        "ibn_label1": "summary of NBI",
        "ibn_label2": "full data of NBI"
    },
    "Português": {
        "title": "Análise DRIS",
        "upload_label": "Carregue seu arquivo de dados (xlsx ou csv)",
        "show_example_label": "Mostrar Arquivo de Exemplo",
        "help_label": "Ajuda",
        "response_variable_label": "Selecione a variável resposta (produtividade):",
        "download_results_label": "Baixar Resultados",
        "dris_indices_label": "Índices DRIS",
        "ibn_label": "IBN",
        "means_and_sds_label": "Médias e Desvios Padrão para Alta Produtividade",
        "dris_indices_high_yield_label": "Índices DRIS para Alta Produtividade",
        "reference_values_label": "Valores de Referência",
        "bootstrap_iterations_label": "Número de iterações Bootstrap",
        "run_bootstrap_label": "Executar Bootstrap",
        "example_file_not_found_error": "O arquivo de exemplo '{example_file}' não foi encontrado.",
        "faqs_title": "Perguntas Frequentes (FAQs)",
        "faqs_q1": "Quais tipos de arquivos posso carregar?",
        "faqs_a1": "Você pode carregar arquivos Excel (.xlsx) ou CSV (.csv) contendo dados de nutrientes e a variável resposta.",
        "faqs_q2": "Como seleciono a variável resposta?",
        "faqs_a2": "Após carregar o arquivo, escolha a coluna que representa a variável resposta (produtividade) no menu dropdown 'Selecione a variável resposta (produtividade):'.",
        "faqs_q3": "O que são os índices DRIS?",
        "faqs_a3": "Os índices DRIS medem a relação entre a concentração de diferentes nutrientes na planta, indicando se a planta está bem nutrida ou se há algum desequilíbrio.",
        "faqs_q4": "O que é o IBN?",
        "faqs_a4": "O IBN (Índice de Balanço Nutricional) é uma medida que resume o grau de desequilíbrio nutricional na planta. Um IBN alto indica um desequilíbrio maior.",
        "faqs_q5": "O que é bootstrap?",
        "faqs_a5": "A análise de bootstrap é uma técnica estatística que permite calcular intervalos de confiança para os índices DRIS e IBN.",
        "faqs_q6": "Como interpreto os valores de referência?",
        "faqs_a6": "Os valores de referência indicam a concentração ideal de cada nutriente na planta para alcançar alta produtividade.",
        "example_usage_title": "Exemplo Detalhado de Uso",
        "example_usage_step1": "1. Carregue o arquivo de dados:",
        "example_usage_step1_download": "Baixe o arquivo de exemplo [link para o arquivo de exemplo].",
        "example_usage_step1_upload": "Carregue o arquivo no aplicativo usando o botão 'Carregue seu arquivo de dados'.",
        "example_usage_step2": "2. Selecione a variável resposta:",
        "example_usage_step2_select": "Selecione a coluna 'Produtividade' como variável resposta.",
        "example_usage_step3": "3. Execute a análise:",
        "example_usage_step3_run": "Clique no botão 'Executar Análise' para calcular os índices DRIS, IBN e outras informações relevantes.",
        "example_usage_step4": "4. Visualize os resultados:",
        "example_usage_step4_view": "Explore as abas 'Índices DRIS', 'Bootstrap DRIS/IBN', 'Médias e Desvios Padrão', 'Índices DRIS Alta Produtividade' e 'Valores de Referência' para ver os resultados da análise.",
        "example_usage_step5": "5. Baixe os resultados:",
        "example_usage_step5_download": "Use os botões de download para salvar os resultados em arquivos Excel.",
        "select_language_label": "Selecione o idioma",
        "density_label": "Densidade",
        "bootkde_label": "Distribuição Bootstrap de ",
        "dowlkde_label": "Baixar gráfico Bootstrap de ",
        "dowlreg_label": "Baixar gráfico de regressão de ",
        "as_PDF_label": "como PDF",
        "as_SVG_label": "como SVG",
        "dldrisridge_label": "Baixar gráfico de Ridge DRIS como PDF",
        "dldrisridgesvg_label": "Baixar gráfico de Ridge DRIS como SVG",
        "reference_label": "Baixar todos os gráficos de valores de referência como PDF",
        "allreg_label": "Baixar todos os gráficos de regressão como PDF",
        "dwlibnpdf_label": "Baixar gráfico de  IBN como PDF",
        "dwlibnsvg_label": "Baixar gráfico de IBN como SVG",
        "dris_indices_label1": "resumo dos  valores DRIS",
        "dris_indices_label2": "dados completos dos valores DRIS",
        "ibn_label1": "resumo de IBN",
        "ibn_label2": "rados completos de IBN"
    },
    "Español": {
        "title": "Análisis DRIS",
        "upload_label": "Suba su archivo de datos (xlsx o csv)",
        "show_example_label": "Mostrar Archivo de Ejemplo",
        "help_label": "Ayuda",
        "response_variable_label": "Seleccione la variable de respuesta (productividad):",
        "download_results_label": "Descargar Resultados",
        "dris_indices_label": "Índices DRIS",
        "ibn_label": "IBN",
        "means_and_sds_label": "Medias y Desviaciones Estándar para Alta Productividad",
        "dris_indices_high_yield_label": "Índices DRIS para Alta Productividad",
        "reference_values_label": "Valores de Referencia",
        "bootstrap_iterations_label": "Número de iteraciones Bootstrap",
        "run_bootstrap_label": "Ejecutar Bootstrap",
        "example_file_not_found_error": "El archivo de ejemplo '{example_file}' no se encontró.",
        "faqs_title": "Preguntas Frecuentes (FAQs)",
        "faqs_q1": "¿Qué tipos de archivos puedo cargar?",
        "faqs_a1": "Puede cargar archivos Excel (.xlsx) o CSV (.csv) que contengan datos de nutrientes y la variable de respuesta.",
        "faqs_q2": "¿Cómo selecciono la variable de respuesta?",
        "faqs_a2": "Después de cargar el archivo, elija la columna que representa la variable de respuesta (productividad) en el menú desplegable 'Seleccione la variable de respuesta (productividad):'.",
        "faqs_q3": "¿Qué son los índices DRIS?",
        "faqs_a3": "Los índices DRIS miden la relación entre la concentración de diferentes nutrientes en la planta, indicando si la planta está bien nutrida o si hay algún desequilibrio.",
        "faqs_q4": "¿Qué es el IBN?",
        "faqs_a4": "El IBN (Índice de Balance Nutricional) es una medida que resume el grado de desequilibrio nutricional en la planta. Un IBN alto indica un desequilibrio mayor.",
        "faqs_q5": "¿Qué es bootstrap?",
        "faqs_a5": "El análisis de bootstrap es una técnica estadística que permite calcular intervalos de confianza para los índices DRIS e IBN.",
        "faqs_q6": "¿Cómo interpreto los valores de referencia?",
        "faqs_a6": "Los valores de referencia indican la concentración ideal de cada nutriente en la planta para alcanzar alta productividad.",
        "example_usage_title": "Ejemplo Detallado de Uso",
        "example_usage_step1": "1. Cargue el archivo de datos:",
        "example_usage_step1_download": "Descargue el archivo de ejemplo [link for example file].",
        "example_usage_step1_upload": "Cargue el archivo en la aplicación usando el botón 'Suba su archivo de datos'.",
        "example_usage_step2": "2. Seleccione la variable de respuesta:",
        "example_usage_step2_select": "Seleccione la columna 'Productividad' como variable de respuesta.",
        "example_usage_step3": "3. Ejecute el análisis:",
        "example_usage_step3_run": "Haga clic en el botón 'Ejecutar Análisis' para calcular los índices DRIS, IBN y otra información relevante.",
        "example_usage_step4": "4. Visualice los resultados:",
        "example_usage_step4_view": "Explore las pestañas 'Índices DRIS', 'Bootstrap DRIS/IBN', 'Medias y Desviaciones Estándar', 'Índices DRIS Alta Productividad' y 'Valores de Referencia' para ver los resultados del análisis.",
        "example_usage_step5": "5. Descargue los resultados:",
        "example_usage_step5_download": "Use los botones de descarga para guardar los resultados en archivos Excel.",
        "select_language_label": "Seleccione el idioma",
        "density_label": "Densidad",
        "bootkde_label": "Distribución Bootstrap de ",
        "dowlkde_label": "Descargar gráfico Bootstrap de ",
        "dowlreg_label": "Descargar gráfico de regresión de ",
        "as_PDF_label": "como PDF",
        "as_SVG_label": "como SVG",
        "dldrisridge_label": "Descargar gráfico de Ridge DRIS como PDF",
        "dldrisridgesvg_label": "Descargar gráfico de Ridge DRIS como SVG",
        "reference_label": "Descargar todos los gráficos de valores de referencia como PDF",
        "allreg_label": "Descargar todos los gráficos de regresión como PDF",
        "dwlibnpdf_label": "Descargar gráfico de  IBN como PDF",
        "dwlibnsvg_label": "Descargar gráfico de IBN como SVG",
        "dris_indices_label1": "resumen de los  valores DRIS",
        "dris_indices_label2": "datos completos de los valores DRIS",
        "ibn_label1": "resumen de IBN",
        "ibn_label2": "datos completos de IBN"
    }
}

st.markdown(
    """
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/tailwindcss@2.2.19/dist/tailwind.min.css">
    """,
    unsafe_allow_html=True
)

def translate(text):
    return translations[language].get(text, text)

uploaded_file = st.file_uploader(translate("upload_label"), type=['xlsx', 'csv'])
example_file = "data/example_data.xlsx"
if st.button(translate("show_example_label")):
    try:
        df_example = pd.read_excel(example_file)
        st.write("## Arquivo de Exemplo:")
        st.dataframe(df_example)
    except FileNotFoundError:
        st.error(translate("example_file_not_found_error").format(example_file=example_file))

if uploaded_file is not None:
    try:
        if uploaded_file.name.endswith('.xlsx'):
            df = pd.read_excel(uploaded_file)
        else:
            df = pd.read_csv(uploaded_file)
        response_variable = st.selectbox(translate("response_variable_label"), df.columns)
        nutrients = [col for col in df.columns if col != response_variable]
        progress_bar = st.progress(0)
        progress_value = 0
        dris_indices, IBN, df_ratios, df_high_yield, all_variance_ratio_ratios = calculate_dris_ibn(df, nutrients, response_variable)
        tab1, tab2, tab3, tab4, tab5 = st.tabs([translate("dris_indices_label"), 
                                                translate("run_bootstrap_label"), 
                                                translate("means_and_sds_label"),
                                                translate("dris_indices_high_yield_label"),
                                                translate("reference_values_label")])
        with tab1:
            st.write(f"## {translate('dris_indices_label')}")
            st.write(pd.DataFrame(dris_indices, index=['DRIS Index']).transpose())
            st.write(f"## {translate('ibn_label')}: {IBN:.2f}")
            st.write(f"**{translate('download_results_label')} {translate('dris_indices_label')}:**")
            download_xlsx(pd.DataFrame(dris_indices, index=['DRIS Index']).transpose(), 'dris_indices.xlsx')
            st.write(f"**{translate('download_results_label')} {translate('ibn_label')}:**") 
            download_xlsx(pd.DataFrame({'IBN': [IBN]}), 'ibn.xlsx')

        
        with tab2:
            st.write(f"## {translate('run_bootstrap_label')}")
            bootstrap_iterations = st.number_input(translate("bootstrap_iterations_label"), min_value=100, value=1000, step=100, key="bootstrap_iterations")
            
            if "tab2_results" not in st.session_state:
                st.session_state.tab2_results = None

            if st.button(translate("run_bootstrap_label"), key="run_bootstrap_tab2"):
                results = process_tab2(df, nutrients, response_variable, bootstrap_iterations)
                st.session_state.tab2_results = results

            if st.session_state.tab2_results:
                results = st.session_state.tab2_results
                dris_indices_samples = [result[0] for result in results]
                IBN_samples = [result[1] for result in results]
                df_dris_indices_samples = pd.DataFrame(dris_indices_samples, columns=nutrients)
                st.write(f"### {translate('run_bootstrap_label')} {translate('dris_indices_label')}:")
                st.write(df_dris_indices_samples.describe(percentiles=[0.025, 0.975]))
                st.write(f"**{translate('download_results_label')} {translate('run_bootstrap_label')} ({translate('dris_indices_label1')}):**")
                download_xlsx(df_dris_indices_samples.describe(percentiles=[0.025, 0.975]), 'bootstrap_dris.xlsx')
                st.write(f"**{translate('download_results_label')} {translate('run_bootstrap_label')} ({translate('dris_indices_label2')}):**")
                download_xlsx(df_dris_indices_samples, 'dris_indices_samples.xlsx')
                st.write(f"### {translate('run_bootstrap_label')} {translate('ibn_label')}:")
                st.write(pd.Series(IBN_samples).describe(percentiles=[0.025, 0.975]))
                st.write(f"**{translate('download_results_label')} {translate('run_bootstrap_label')} ({translate('ibn_label1')}):**") 
                download_xlsx(pd.DataFrame(pd.Series(IBN_samples).describe(percentiles=[0.025, 0.975])), 'bootstrap_ibn.xlsx')
                st.write(f"**{translate('download_results_label')} {translate('run_bootstrap_label')} ({translate('ibn_label2')}):**")
                download_xlsx(pd.DataFrame({'IBN_samples': IBN_samples}), 'ibn_samples.xlsx')
                sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
                # Derreter o DataFrame para trabalhar com o seaborn FacetGrid
                # Derreter o DataFrame para trabalhar com o seaborn FacetGrid
                df_dris_melted = df_dris_indices_samples.melt(var_name='nutrient', value_name='dris_index')

                # Encontrar os limites globais para o eixo x
                global_lower_bound = df_dris_indices_samples.quantile(0.025).min()
                global_upper_bound = df_dris_indices_samples.quantile(0.975).max()

                # Criar o FacetGrid
                g = sns.FacetGrid(df_dris_melted, row='nutrient', hue='nutrient', aspect=15, height=.5, palette='viridis')

                # Mapear o KDE plot para o grid
                g.map(sns.kdeplot, 'dris_index', bw_adjust=.5, clip_on=False, fill=True, alpha=0.5, linewidth=1.5, color='green')

                # Mapear o KDE plot para criar a linha branca
                g.map(sns.kdeplot, 'dris_index', clip_on=False, color="w", lw=2, bw_adjust=.5)

                # Adicionar uma linha de referência
                g.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)

                # Função para adicionar labels
                def label(x, color, label):
                    ax = plt.gca()
                    ax.text(0, .2, label, fontweight="bold", color='black',
                            ha="left", va="center", transform=ax.transAxes)

                # Mapear a função label para adicionar os labels
                g.map(label, 'dris_index')

                # Ajustar o espaço entre os subplots
                g.figure.subplots_adjust(hspace=-.25)

                # Definir os títulos e labels dos eixos
                g.set_titles("")
                g.set(yticks=[], ylabel="")

                # Remover as bordas desnecessárias
                g.despine(bottom=True, left=True)

                # Ajustar o eixo x de cada subplot individualmente com linhas de intervalo de confiança
                for i, nutrient in enumerate(df_dris_indices_samples.columns):
                    lower_bound = df_dris_indices_samples[nutrient].quantile(0.025)
                    upper_bound = df_dris_indices_samples[nutrient].quantile(0.975)
                    
                    # Adicionar linhas de intervalo de confiança individuais
                    g.axes.flat[i].axvline(x=lower_bound, color='red', linestyle='--', linewidth=1)
                    g.axes.flat[i].axvline(x=upper_bound, color='red', linestyle='--', linewidth=1)
                    
                    # Setar o limite global do eixo x
                    g.axes.flat[i].set_xlim(global_lower_bound, global_upper_bound)
                    
                # Definir o título do eixo x e y
                g.set_axis_labels("DRIS", "")

                # Mostrar o gráfico no Streamlit
                st.pyplot(g.fig)
                pdf_buffer = io.BytesIO()
                with PdfPages(pdf_buffer) as pdf:
                    pdf.savefig(g.fig)
                st.download_button(label=translate("dldrisridge_label"),
                                data=pdf_buffer.getvalue(),
                                file_name="dris_bootstrap_ridge_plots.pdf",
                                mime="application/pdf")
                svg_buffer = io.BytesIO()
                g.fig.savefig(svg_buffer, format='svg')
                st.download_button(label=translate("dldrisridgesvg_label"),
                                data=svg_buffer.getvalue(),
                                file_name="dris_bootstrap_ridge_plots.svg",
                                mime="image/svg+xml")
                plt.close(g.fig)
                plt.figure(figsize=(8, 4))
                sns.kdeplot(IBN_samples, color='green', fill=True)
                plt.xlabel(translate("ibn_label"))
                plt.ylabel(translate("density_label"))
                plt.title(translate("bootkde_label") + translate("ibn_label"))
                lower_bound_ibn = np.percentile(IBN_samples, 2.5)
                upper_bound_ibn = np.percentile(IBN_samples, 97.5)
                plt.axvline(lower_bound_ibn, color='red', linestyle='--', linewidth=1)
                plt.axvline(upper_bound_ibn, color='red', linestyle='--', linewidth=1)
                st.pyplot(plt.gcf())
                ibn_pdf_buffer = io.BytesIO()
                with PdfPages(ibn_pdf_buffer) as pdf:
                    pdf.savefig(plt.gcf())
                st.download_button(label=translate("dwlibnpdf_label"),
                                data=ibn_pdf_buffer.getvalue(),
                                file_name="ibn_bootstrap_plot.pdf",
                                mime="application/pdf")
                ibn_svg_buffer = io.BytesIO()
                plt.savefig(ibn_svg_buffer, format='svg')
                st.download_button(label=translate("dwlibnsvg_label"),
                                data=ibn_svg_buffer.getvalue(),
                                file_name="ibn_bootstrap_plot.svg",
                                mime="image/svg+xml")
                plt.close()
                plt.savefig('dris_bootstrap_plot_intervalos.png', dpi=600)
                plt.close()
                progress_value += 20
                progress_bar.progress(progress_value)
        with tab3:
            st.write(f"## {translate('means_and_sds_label')}")
            st.write(df_high_yield[nutrients].describe(percentiles=[0.025, 0.975]))
            st.write(f"**{translate('download_results_label')} ({translate('means_and_sds_label')}):**")
            download_xlsx(df_high_yield[nutrients].describe(percentiles=[0.025, 0.975]), 'high_yield_stats.xlsx')
            st.write(f"**{translate('download_results_label')} ({translate('means_and_sds_label')}):**")
            download_xlsx(df_high_yield, 'df_high_yield.xlsx')
        with tab4:
            st.write(f"## {translate('dris_indices_high_yield_label')}")
            dris_indices_high_yield = {}
            for nutrient in nutrients:
                f_high_total = 0
                f_low_total = 0
                for ratio, variance_ratio in all_variance_ratio_ratios[nutrient].items():
                    numer, denom = ratio.split('/')
                    if numer == nutrient:
                        Ea_num = (df_high_yield[numer] / df_high_yield[denom])
                        En_num = (df_high_yield[numer] / df_high_yield[denom]).mean()
                        s_num = (df_high_yield[numer] / df_high_yield[denom]).std()
                        f_high = f_E_Ei(Ea_num, En_num, s_num)
                        f_high_total += f_high
                    else:
                        Ea_den = (df_high_yield[numer] / df_high_yield[nutrient])
                        En_den = (df_high_yield[numer] / df_high_yield[nutrient]).mean()
                        s_den = (df_high_yield[numer] / df_high_yield[nutrient]).std()
                        f_low = f_E_Ei(Ea_den, En_den, s_den)
                        f_low_total += f_low
                dris_index = (f_high_total - f_low_total) / (len(nutrients) - 1)
                dris_indices_high_yield[nutrient] = dris_indices_high_yield.get(nutrient, 0) + dris_index
            dris_indices_high_yield_df = pd.DataFrame(dris_indices_high_yield)
            st.write(dris_indices_high_yield_df)
            download_xlsx(dris_indices_high_yield_df, 'dris_indices_high_yield.xlsx')
        with tab5:
            zeros = {}
            regression_models = {}
            for nutrient in nutrients:
                y = df_high_yield[[nutrient]].values.flatten()
                X = dris_indices_high_yield_df[[nutrient]].values.flatten()
                coefs = np.polyfit(X, y, 1)
                intercept = coefs[1]
                slope = coefs[0]
                regression_models[nutrient] = (slope, intercept)
                zero = intercept
                std = y.std()
                interval = [zero - (2 / 3) * std, zero + (2 / 3) * std]
                zeros[nutrient] = (zero, interval)
            model_parameters_df = pd.DataFrame(
                regression_models, index=["Slope", "Intercept"]
            ).T
            model_parameters_df["R^2"] = [
                np.corrcoef(dris_indices_high_yield_df[nutrient], df_high_yield[nutrient])[0, 1] ** 2
                for nutrient in nutrients
            ]
            st.write(model_parameters_df)
            download_xlsx(model_parameters_df, 'model_parameters.xlsx')
            bootstrap_iterations = st.number_input(
                translate("bootstrap_iterations_label"),
                min_value=100,
                value=1000,
                step=100,
                key="bootstrap_iterations_reference"
            )
            if "tab5_results" not in st.session_state:
                st.session_state.tab5_results = None

            if st.button(translate("run_bootstrap_label"), key="run_bootstrap_tab5"):
                df_zeros_samples = process_tab5(df_high_yield, dris_indices_high_yield_df, nutrients, bootstrap_iterations)
                st.session_state.tab5_results = df_zeros_samples

            if st.session_state.tab5_results is not None:
                df_zeros_samples = st.session_state.tab5_results
                st.write(df_zeros_samples.describe(percentiles=[0.025, 0.975]))
                download_xlsx(df_zeros_samples.describe(percentiles=[0.025, 0.975]), 'bootstrap_reference_values.xlsx')
                download_xlsx(df_zeros_samples, 'df_zeros_samples.xlsx')
                
                all_plots_pdf_buffer = io.BytesIO()
                with PdfPages(all_plots_pdf_buffer) as all_plots_pdf:
                    for nutrient in nutrients:
                        plot_kde(nutrient, df_zeros_samples, pdf=all_plots_pdf)

                st.download_button(label=translate("reference_label"), key="dwltab5_all_pdf",
                                data=all_plots_pdf_buffer.getvalue(),
                                file_name="all_bootstrap_ridge_plots.pdf",
                                mime="application/pdf")
                
                regression_pdf_buffer = io.BytesIO()
                regression_pdf = PdfPages(regression_pdf_buffer)
                for nutrient in nutrients:
                    plot_original_regression(nutrient)
                regression_pdf.close()
                st.download_button(label=translate("allreg_label"), key="dwltab52_original",
                                data=regression_pdf_buffer.getvalue(),
                                file_name="original_regression_plots.pdf",
                                mime="application/pdf")
        progress_value += 20
        progress_bar.progress(progress_value)
    except Exception as e:
        st.error(f"Erro ao carregar o arquivo: {e}")

if st.button(translate("help_label")):
    st.markdown(
        f"""
        # {translate("title")} - Streamlit App
        {translate("example_usage_title")}
        {translate("example_usage_step1")}
        {translate("example_usage_step1_download")}
        {translate("example_usage_step1_upload")}
        {translate("example_usage_step2")}
        {translate("example_usage_step2_select")}
        {translate("example_usage_step3")}
        {translate("example_usage_step3_run")}
        {translate("example_usage_step4")}
        {translate("example_usage_step4_view")}
        {translate("example_usage_step5")}
        {translate("example_usage_step5_download")}
        ## {translate("faqs_title")}
        **1. {translate("faqs_q1")}**
        {translate("faqs_a1")}
        **2. {translate("faqs_q2")}**
        {translate("faqs_a2")}
        **3. {translate("faqs_q3")}**
        {translate("faqs_a3")}
        **4. {translate("faqs_q4")}**
        {translate("faqs_a4")}
        **5. {translate("faqs_q5")}**
        {translate("faqs_a5")}
        **6. {translate("faqs_q6")}**
        {translate("faqs_a6")}
        """,
        unsafe_allow_html=True
    )
