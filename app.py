import streamlit as st
import pandas as pd
from itertools import permutations
import numpy as np
import io


# Função para download de DataFrame como xlsx
def download_xlsx(df, filename):
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name='Sheet1', index=False)
    writer.close()  
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

    # Calcula todas as razões primeiro
    all_ratios = {}
    for pair in pairs:
        ratio_name = f'{pair[0]}/{pair[1]}'
        all_ratios[ratio_name] = df[pair[0]] / df[pair[1]]

    # Une as razões ao DataFrame df_ratios
    df_ratios = pd.concat([df, pd.DataFrame(all_ratios)], axis=1)

    # Calcular o limite para alto rendimento
    threshold = df[response_variable].mean() + 0.5 * df[response_variable].std()

    # Criar data frames para populações de alto e baixo rendimento
    df_high_yield = df[df[response_variable] >= threshold]
    df_low_yield = df[df[response_variable] < threshold]

    # Calcular as razões de nutrientes para data frames de alto e baixo rendimento separadamente
    df_high_yield_ratios = df_high_yield.copy()
    df_low_yield_ratios = df_low_yield.copy()
    for pair in pairs:
        df_high_yield_ratios[f"{pair[0]}/{pair[1]}"] = df_high_yield_ratios[pair[0]] / df_high_yield_ratios[pair[1]]
        df_low_yield_ratios[f"{pair[0]}/{pair[1]}"] = df_low_yield_ratios[pair[0]] / df_low_yield_ratios[pair[1]]

    # Calcular a variância para cada razão de nutriente em ambas as populações
    var_high_yield_ratios = df_high_yield_ratios.iloc[:, -len(pairs):].var()
    var_low_yield_ratios = df_low_yield_ratios.iloc[:, -len(pairs):].var()

    # Calcular a razão das variâncias
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

# Função para o cálculo de f(E/Ei)
def f_E_Ei(E_Ea, E_En, s, K=10):
    return (E_Ea - E_En) * (K / s)

# Interface do Streamlit
st.title("Análise DRIS")
# Informações da versão e autor
st.sidebar.markdown("**Versão:** 0.1")
st.sidebar.markdown("**Autor:** Prof. Walter E. Pereira, UFPB, CCA")

# Upload de dados
uploaded_file = st.file_uploader("Carregue seu arquivo de dados (xlsx ou csv)", type=['xlsx', 'csv'])
if uploaded_file is not None:
    try:
        if uploaded_file.name.endswith('.xlsx'):
            df = pd.read_excel(uploaded_file)
        else:  # Assume CSV
            df = pd.read_csv(uploaded_file)
        
        # Selecionar a variável resposta (produtividade)
        response_variable = st.selectbox("Selecione a variável resposta (produtividade):", df.columns)

        # Lista de nutrientes (todas as colunas menos a variável resposta)
        nutrients = [col for col in df.columns if col != response_variable]

        # Calcular os índices DRIS e IBN
        dris_indices, IBN, df_ratios, df_high_yield, all_variance_ratio_ratios = calculate_dris_ibn(df, nutrients, response_variable)

        # Criar as abas
        tab1, tab2, tab3, tab4, tab5 = st.tabs(["Índices DRIS", "Bootstrap DRIS/IBN", "Médias e Desvios Padrão",
                                                "Índices DRIS Alta Produtividade", "Valores de Referência"])

        # Aba 1: Índices DRIS
        with tab1:
            st.write("## Índices DRIS")
            st.write(pd.DataFrame(dris_indices, index=['DRIS Index']).transpose())
            st.write(f"## IBN: {IBN:.2f}")
            
            st.write("**Download dos Índices DRIS:**")
            download_xlsx(pd.DataFrame(dris_indices, index=['DRIS Index']).transpose(), 'dris_indices.xlsx')
            
            st.write("**Download do IBN:**") 
            download_xlsx(pd.DataFrame({'IBN': [IBN]}), 'ibn.xlsx') 

        # Aba 2: Bootstrap DRIS/IBN
        with tab2:
            st.write("## Bootstrap DRIS/IBN")
            bootstrap_iterations = st.number_input("Número de iterações Bootstrap:", min_value=100, value=1000, step=100)
            if st.button("Executar Bootstrap"):
                dris_indices_samples = []
                IBN_samples = []
                for _ in range(bootstrap_iterations):
                    df_sample = df.sample(n=len(df), replace=True)
                    dris_indices_sample, IBN_sample, _, _, _ = calculate_dris_ibn(df_sample, nutrients, response_variable)
                    dris_indices_samples.append(dris_indices_sample)
                    IBN_samples.append(IBN_sample)

                df_dris_indices_samples = pd.DataFrame(dris_indices_samples, columns=nutrients)
                st.write("### Resultados do Bootstrap para Índices DRIS:")
                st.write(df_dris_indices_samples.describe(percentiles=[0.025, 0.975]))
                st.write("**Download dos Resultados do Bootstrap (Índices DRIS):**")
                download_xlsx(df_dris_indices_samples.describe(percentiles=[0.025, 0.975]), 'bootstrap_dris.xlsx')
                st.write("**Download das Amostras do Bootstrap (Índices DRIS):**")
                download_xlsx(df_dris_indices_samples, 'dris_indices_samples.xlsx')

                st.write("### Resultados do Bootstrap para IBN:")
                st.write(pd.Series(IBN_samples).describe(percentiles=[0.025, 0.975]))
                st.write("**Download dos Resultados do Bootstrap (IBN):**") 
                download_xlsx(pd.DataFrame(pd.Series(IBN_samples).describe(percentiles=[0.025, 0.975])), 'bootstrap_ibn.xlsx')
                st.write("**Download das Amostras do Bootstrap (IBN):**")
                download_xlsx(pd.DataFrame({'IBN_samples': IBN_samples}), 'ibn_samples.xlsx') 

        # Aba 3: Médias e Desvios Padrão
        with tab3:
            st.write("## Médias e Desvios Padrão para Alta Produtividade")
            st.write(df_high_yield[nutrients].describe(percentiles=[0.025, 0.975]))
            st.write("**Download das Estatísticas Descritivas (Alta Produtividade):**")
            download_xlsx(df_high_yield[nutrients].describe(percentiles=[0.025, 0.975]), 'high_yield_stats.xlsx')
            st.write("**Download do DataFrame de Alta Produtividade:**")
            download_xlsx(df_high_yield, 'df_high_yield.xlsx') 

        # Aba 4: Índices DRIS Alta Produtividade
        with tab4:
            st.write("## Índices DRIS para Alta Produtividade")
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
        
        # Aba 5: Valores de Referência
        with tab5:
            zeros = {}
            regression_models = {}

            for nutrient in nutrients:
                # Calcular a regressão usando np.polyfit
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

            # Bootstrap para intervalos de confiança dos Valores de Referência
            bootstrap_iterations = st.number_input(
                "Número de iterações Bootstrap (Valores de Referência):",
                min_value=100,
                value=1000,
                step=100,
            )

            if st.button("Executar Bootstrap (Valores de Referência)"):
                zeros_samples = []
                for _ in range(bootstrap_iterations):
                    zeros_sample = {}
                    for nutrient in nutrients:
                        df_sample = df_high_yield.sample(n=len(df_high_yield), replace=True)
                        dris_sample = pd.DataFrame({nutrient: dris_indices_high_yield_df[nutrient].sample(n=len(dris_indices_high_yield_df), replace=True)})
                        coeffs = np.polyfit(dris_sample[nutrient], df_sample[nutrient], 1)
                        zero = coeffs[1]
                        zeros_sample[nutrient] = zero
                    zeros_samples.append(zeros_sample)
                df_zeros_samples = pd.DataFrame(zeros_samples)

                st.write(df_zeros_samples.describe(percentiles=[0.025, 0.975]))    
                st.write("**Download dos Valores de Referência do Bootstrap:**")
                download_xlsx(df_zeros_samples.describe(percentiles=[0.025, 0.975]), 'bootstrap_reference_values.xlsx')
                st.write("**Download das Amostras do Bootstrap (Valores de Referência):**")
                download_xlsx(df_zeros_samples, 'df_zeros_samples.xlsx')

    except Exception as e:
        st.error(f"Erro ao carregar o arquivo: {e}")
