{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "id": "de8e274d",
   "metadata": {
    "_cell_guid": "b1076dfc-b9ad-4769-8c92-a6c4dae69d19",
    "_uuid": "8f2839f25d086af736a60e9eeb907d3b93b6e0e5",
    "execution": {
     "iopub.execute_input": "2022-07-30T02:46:17.175012Z",
     "iopub.status.busy": "2022-07-30T02:46:17.172097Z",
     "iopub.status.idle": "2022-07-30T03:00:27.808146Z",
     "shell.execute_reply": "2022-07-30T03:00:27.805851Z"
    },
    "papermill": {
     "duration": 850.645475,
     "end_time": "2022-07-30T03:00:27.810993",
     "exception": false,
     "start_time": "2022-07-30T02:46:17.165518",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Installing package into ‘/usr/local/lib/R/site-library’\n",
      "(as ‘lib’ is unspecified)\n",
      "\n"
     ]
    }
   ],
   "source": [
    "install.packages(\"rstan\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "id": "a1a91b62",
   "metadata": {
    "collapsed": true,
    "execution": {
     "iopub.execute_input": "2022-07-30T03:00:27.950472Z",
     "iopub.status.busy": "2022-07-30T03:00:27.875096Z",
     "iopub.status.idle": "2022-07-30T03:01:17.941656Z",
     "shell.execute_reply": "2022-07-30T03:01:17.939382Z"
    },
    "jupyter": {
     "outputs_hidden": true
    },
    "papermill": {
     "duration": 50.131676,
     "end_time": "2022-07-30T03:01:17.946255",
     "exception": false,
     "start_time": "2022-07-30T03:00:27.814579",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Loading required package: StanHeaders\n",
      "\n",
      "Loading required package: ggplot2\n",
      "\n",
      "rstan (Version 2.21.5, GitRev: 2e1f913d3ca3)\n",
      "\n",
      "For execution on a local, multicore CPU with excess RAM we recommend calling\n",
      "options(mc.cores = parallel::detectCores()).\n",
      "To avoid recompilation of unchanged Stan programs, we recommend calling\n",
      "rstan_options(auto_write = TRUE)\n",
      "\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "\n",
      "stn_md> stancode <- 'data {real y_mean;} parameters {real y;} model {y ~ normal(y_mean,1);}'\n",
      "\n",
      "stn_md> mod <- stan_model(model_code = stancode, verbose = TRUE)\n",
      "\n",
      "TRANSLATING MODEL '16a540c6086086816528e4524def24d9' FROM Stan CODE TO C++ CODE NOW.\n",
      "successful in parsing the Stan model '16a540c6086086816528e4524def24d9'.\n",
      "OS: x86_64, linux-gnu; rstan: 2.21.5; Rcpp: 1.0.9; inline: 0.3.19 \n",
      " >> setting environment variables: \n",
      "PKG_LIBS =  '/usr/local/lib/R/site-library/rstan/lib//libStanServices.a' -L'/usr/local/lib/R/site-library/StanHeaders/lib/' -lStanHeaders -L'/usr/local/lib/R/site-library/RcppParallel/lib/' -ltbb\n",
      "PKG_CPPFLAGS =   -I\"/usr/local/lib/R/site-library/Rcpp/include/\"  -I\"/usr/local/lib/R/site-library/RcppEigen/include/\"  -I\"/usr/local/lib/R/site-library/RcppEigen/include/unsupported\"  -I\"/usr/local/lib/R/site-library/BH/include\" -I\"/usr/local/lib/R/site-library/StanHeaders/include/src/\"  -I\"/usr/local/lib/R/site-library/StanHeaders/include/\"  -I\"/usr/local/lib/R/site-library/RcppParallel/include/\"  -I\"/usr/local/lib/R/site-library/rstan/include\" -DEIGEN_NO_DEBUG  -DBOOST_DISABLE_ASSERTS  -DBOOST_PENDING_INTEGER_LOG2_HPP  -DSTAN_THREADS  -DBOOST_NO_AUTO_PTR  -include '/usr/local/lib/R/site-library/StanHeaders/include/stan/math/prim/mat/fun/Eigen.hpp'  -D_REENTRANT -DRCPP_PARALLEL_USE_TBB=1 \n",
      " >> Program source :\n",
      "\n",
      "   1 : \n",
      "   2 : // includes from the plugin\n",
      "   3 : // [[Rcpp::plugins(cpp14)]]\n",
      "   4 : \n",
      "   5 : \n",
      "   6 : // user includes\n",
      "   7 : #include <Rcpp.h>\n",
      "   8 : #include <rstan/io/rlist_ref_var_context.hpp>\n",
      "   9 : #include <rstan/io/r_ostream.hpp>\n",
      "  10 : #include <rstan/stan_args.hpp>\n",
      "  11 : #include <boost/integer/integer_log2.hpp>\n",
      "  12 : // Code generated by Stan version 2.21.0\n",
      "  13 : \n",
      "  14 : #include <stan/model/model_header.hpp>\n",
      "  15 : \n",
      "  16 : namespace modele188e6315_16a540c6086086816528e4524def24d9_namespace {\n",
      "  17 : \n",
      "  18 : using std::istream;\n",
      "  19 : using std::string;\n",
      "  20 : using std::stringstream;\n",
      "  21 : using std::vector;\n",
      "  22 : using stan::io::dump;\n",
      "  23 : using stan::math::lgamma;\n",
      "  24 : using stan::model::prob_grad;\n",
      "  25 : using namespace stan::math;\n",
      "  26 : \n",
      "  27 : static int current_statement_begin__;\n",
      "  28 : \n",
      "  29 : stan::io::program_reader prog_reader__() {\n",
      "  30 :     stan::io::program_reader reader;\n",
      "  31 :     reader.add_event(0, 0, \"start\", \"modele188e6315_16a540c6086086816528e4524def24d9\");\n",
      "  32 :     reader.add_event(3, 1, \"end\", \"modele188e6315_16a540c6086086816528e4524def24d9\");\n",
      "  33 :     return reader;\n",
      "  34 : }\n",
      "  35 : \n",
      "  36 : class modele188e6315_16a540c6086086816528e4524def24d9\n",
      "  37 :   : public stan::model::model_base_crtp<modele188e6315_16a540c6086086816528e4524def24d9> {\n",
      "  38 : private:\n",
      "  39 :         double y_mean;\n",
      "  40 : public:\n",
      "  41 :     modele188e6315_16a540c6086086816528e4524def24d9(rstan::io::rlist_ref_var_context& context__,\n",
      "  42 :         std::ostream* pstream__ = 0)\n",
      "  43 :         : model_base_crtp(0) {\n",
      "  44 :         ctor_body(context__, 0, pstream__);\n",
      "  45 :     }\n",
      "  46 : \n",
      "  47 :     modele188e6315_16a540c6086086816528e4524def24d9(stan::io::var_context& context__,\n",
      "  48 :         unsigned int random_seed__,\n",
      "  49 :         std::ostream* pstream__ = 0)\n",
      "  50 :         : model_base_crtp(0) {\n",
      "  51 :         ctor_body(context__, random_seed__, pstream__);\n",
      "  52 :     }\n",
      "  53 : \n",
      "  54 :     void ctor_body(stan::io::var_context& context__,\n",
      "  55 :                    unsigned int random_seed__,\n",
      "  56 :                    std::ostream* pstream__) {\n",
      "  57 :         typedef double local_scalar_t__;\n",
      "  58 : \n",
      "  59 :         boost::ecuyer1988 base_rng__ =\n",
      "  60 :           stan::services::util::create_rng(random_seed__, 0);\n",
      "  61 :         (void) base_rng__;  // suppress unused var warning\n",
      "  62 : \n",
      "  63 :         current_statement_begin__ = -1;\n",
      "  64 : \n",
      "  65 :         static const char* function__ = \"modele188e6315_16a540c6086086816528e4524def24d9_namespace::modele188e6315_16a540c6086086816528e4524def24d9\";\n",
      "  66 :         (void) function__;  // dummy to suppress unused var warning\n",
      "  67 :         size_t pos__;\n",
      "  68 :         (void) pos__;  // dummy to suppress unused var warning\n",
      "  69 :         std::vector<int> vals_i__;\n",
      "  70 :         std::vector<double> vals_r__;\n",
      "  71 :         local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());\n",
      "  72 :         (void) DUMMY_VAR__;  // suppress unused var warning\n",
      "  73 : \n",
      "  74 :         try {\n",
      "  75 :             // initialize data block variables from context__\n",
      "  76 :             current_statement_begin__ = 1;\n",
      "  77 :             context__.validate_dims(\"data initialization\", \"y_mean\", \"double\", context__.to_vec());\n",
      "  78 :             y_mean = double(0);\n",
      "  79 :             vals_r__ = context__.vals_r(\"y_mean\");\n",
      "  80 :             pos__ = 0;\n",
      "  81 :             y_mean = vals_r__[pos__++];\n",
      "  82 : \n",
      "  83 : \n",
      "  84 :             // initialize transformed data variables\n",
      "  85 :             // execute transformed data statements\n",
      "  86 : \n",
      "  87 :             // validate transformed data\n",
      "  88 : \n",
      "  89 :             // validate, set parameter ranges\n",
      "  90 :             num_params_r__ = 0U;\n",
      "  91 :             param_ranges_i__.clear();\n",
      "  92 :             current_statement_begin__ = 1;\n",
      "  93 :             num_params_r__ += 1;\n",
      "  94 :         } catch (const std::exception& e) {\n",
      "  95 :             stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());\n",
      "  96 :             // Next line prevents compiler griping about no return\n",
      "  97 :             throw std::runtime_error(\"*** IF YOU SEE THIS, PLEASE REPORT A BUG ***\");\n",
      "  98 :         }\n",
      "  99 :     }\n",
      " 100 : \n",
      " 101 :     ~modele188e6315_16a540c6086086816528e4524def24d9() { }\n",
      " 102 : \n",
      " 103 : \n",
      " 104 :     void transform_inits(const stan::io::var_context& context__,\n",
      " 105 :                          std::vector<int>& params_i__,\n",
      " 106 :                          std::vector<double>& params_r__,\n",
      " 107 :                          std::ostream* pstream__) const {\n",
      " 108 :         typedef double local_scalar_t__;\n",
      " 109 :         stan::io::writer<double> writer__(params_r__, params_i__);\n",
      " 110 :         size_t pos__;\n",
      " 111 :         (void) pos__; // dummy call to supress warning\n",
      " 112 :         std::vector<double> vals_r__;\n",
      " 113 :         std::vector<int> vals_i__;\n",
      " 114 : \n",
      " 115 :         current_statement_begin__ = 1;\n",
      " 116 :         if (!(context__.contains_r(\"y\")))\n",
      " 117 :             stan::lang::rethrow_located(std::runtime_error(std::string(\"Variable y missing\")), current_statement_begin__, prog_reader__());\n",
      " 118 :         vals_r__ = context__.vals_r(\"y\");\n",
      " 119 :         pos__ = 0U;\n",
      " 120 :         context__.validate_dims(\"parameter initialization\", \"y\", \"double\", context__.to_vec());\n",
      " 121 :         double y(0);\n",
      " 122 :         y = vals_r__[pos__++];\n",
      " 123 :         try {\n",
      " 124 :             writer__.scalar_unconstrain(y);\n",
      " 125 :         } catch (const std::exception& e) {\n",
      " 126 :             stan::lang::rethrow_located(std::runtime_error(std::string(\"Error transforming variable y: \") + e.what()), current_statement_begin__, prog_reader__());\n",
      " 127 :         }\n",
      " 128 : \n",
      " 129 :         params_r__ = writer__.data_r();\n",
      " 130 :         params_i__ = writer__.data_i();\n",
      " 131 :     }\n",
      " 132 : \n",
      " 133 :     void transform_inits(const stan::io::var_context& context,\n",
      " 134 :                          Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,\n",
      " 135 :                          std::ostream* pstream__) const {\n",
      " 136 :       std::vector<double> params_r_vec;\n",
      " 137 :       std::vector<int> params_i_vec;\n",
      " 138 :       transform_inits(context, params_i_vec, params_r_vec, pstream__);\n",
      " 139 :       params_r.resize(params_r_vec.size());\n",
      " 140 :       for (int i = 0; i < params_r.size(); ++i)\n",
      " 141 :         params_r(i) = params_r_vec[i];\n",
      " 142 :     }\n",
      " 143 : \n",
      " 144 : \n",
      " 145 :     template <bool propto__, bool jacobian__, typename T__>\n",
      " 146 :     T__ log_prob(std::vector<T__>& params_r__,\n",
      " 147 :                  std::vector<int>& params_i__,\n",
      " 148 :                  std::ostream* pstream__ = 0) const {\n",
      " 149 : \n",
      " 150 :         typedef T__ local_scalar_t__;\n",
      " 151 : \n",
      " 152 :         local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());\n",
      " 153 :         (void) DUMMY_VAR__;  // dummy to suppress unused var warning\n",
      " 154 : \n",
      " 155 :         T__ lp__(0.0);\n",
      " 156 :         stan::math::accumulator<T__> lp_accum__;\n",
      " 157 :         try {\n",
      " 158 :             stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);\n",
      " 159 : \n",
      " 160 :             // model parameters\n",
      " 161 :             current_statement_begin__ = 1;\n",
      " 162 :             local_scalar_t__ y;\n",
      " 163 :             (void) y;  // dummy to suppress unused var warning\n",
      " 164 :             if (jacobian__)\n",
      " 165 :                 y = in__.scalar_constrain(lp__);\n",
      " 166 :             else\n",
      " 167 :                 y = in__.scalar_constrain();\n",
      " 168 : \n",
      " 169 :             // model body\n",
      " 170 : \n",
      " 171 :             current_statement_begin__ = 1;\n",
      " 172 :             lp_accum__.add(normal_log<propto__>(y, y_mean, 1));\n",
      " 173 : \n",
      " 174 :         } catch (const std::exception& e) {\n",
      " 175 :             stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());\n",
      " 176 :             // Next line prevents compiler griping about no return\n",
      " 177 :             throw std::runtime_error(\"*** IF YOU SEE THIS, PLEASE REPORT A BUG ***\");\n",
      " 178 :         }\n",
      " 179 : \n",
      " 180 :         lp_accum__.add(lp__);\n",
      " 181 :         return lp_accum__.sum();\n",
      " 182 : \n",
      " 183 :     } // log_prob()\n",
      " 184 : \n",
      " 185 :     template <bool propto, bool jacobian, typename T_>\n",
      " 186 :     T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,\n",
      " 187 :                std::ostream* pstream = 0) const {\n",
      " 188 :       std::vector<T_> vec_params_r;\n",
      " 189 :       vec_params_r.reserve(params_r.size());\n",
      " 190 :       for (int i = 0; i < params_r.size(); ++i)\n",
      " 191 :         vec_params_r.push_back(params_r(i));\n",
      " 192 :       std::vector<int> vec_params_i;\n",
      " 193 :       return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);\n",
      " 194 :     }\n",
      " 195 : \n",
      " 196 : \n",
      " 197 :     void get_param_names(std::vector<std::string>& names__) const {\n",
      " 198 :         names__.resize(0);\n",
      " 199 :         names__.push_back(\"y\");\n",
      " 200 :     }\n",
      " 201 : \n",
      " 202 : \n",
      " 203 :     void get_dims(std::vector<std::vector<size_t> >& dimss__) const {\n",
      " 204 :         dimss__.resize(0);\n",
      " 205 :         std::vector<size_t> dims__;\n",
      " 206 :         dims__.resize(0);\n",
      " 207 :         dimss__.push_back(dims__);\n",
      " 208 :     }\n",
      " 209 : \n",
      " 210 :     template <typename RNG>\n",
      " 211 :     void write_array(RNG& base_rng__,\n",
      " 212 :                      std::vector<double>& params_r__,\n",
      " 213 :                      std::vector<int>& params_i__,\n",
      " 214 :                      std::vector<double>& vars__,\n",
      " 215 :                      bool include_tparams__ = true,\n",
      " 216 :                      bool include_gqs__ = true,\n",
      " 217 :                      std::ostream* pstream__ = 0) const {\n",
      " 218 :         typedef double local_scalar_t__;\n",
      " 219 : \n",
      " 220 :         vars__.resize(0);\n",
      " 221 :         stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);\n",
      " 222 :         static const char* function__ = \"modele188e6315_16a540c6086086816528e4524def24d9_namespace::write_array\";\n",
      " 223 :         (void) function__;  // dummy to suppress unused var warning\n",
      " 224 : \n",
      " 225 :         // read-transform, write parameters\n",
      " 226 :         double y = in__.scalar_constrain();\n",
      " 227 :         vars__.push_back(y);\n",
      " 228 : \n",
      " 229 :         double lp__ = 0.0;\n",
      " 230 :         (void) lp__;  // dummy to suppress unused var warning\n",
      " 231 :         stan::math::accumulator<double> lp_accum__;\n",
      " 232 : \n",
      " 233 :         local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());\n",
      " 234 :         (void) DUMMY_VAR__;  // suppress unused var warning\n",
      " 235 : \n",
      " 236 :         if (!include_tparams__ && !include_gqs__) return;\n",
      " 237 : \n",
      " 238 :         try {\n",
      " 239 :             if (!include_gqs__ && !include_tparams__) return;\n",
      " 240 :             if (!include_gqs__) return;\n",
      " 241 :         } catch (const std::exception& e) {\n",
      " 242 :             stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());\n",
      " 243 :             // Next line prevents compiler griping about no return\n",
      " 244 :             throw std::runtime_error(\"*** IF YOU SEE THIS, PLEASE REPORT A BUG ***\");\n",
      " 245 :         }\n",
      " 246 :     }\n",
      " 247 : \n",
      " 248 :     template <typename RNG>\n",
      " 249 :     void write_array(RNG& base_rng,\n",
      " 250 :                      Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,\n",
      " 251 :                      Eigen::Matrix<double,Eigen::Dynamic,1>& vars,\n",
      " 252 :                      bool include_tparams = true,\n",
      " 253 :                      bool include_gqs = true,\n",
      " 254 :                      std::ostream* pstream = 0) const {\n",
      " 255 :       std::vector<double> params_r_vec(params_r.size());\n",
      " 256 :       for (int i = 0; i < params_r.size(); ++i)\n",
      " 257 :         params_r_vec[i] = params_r(i);\n",
      " 258 :       std::vector<double> vars_vec;\n",
      " 259 :       std::vector<int> params_i_vec;\n",
      " 260 :       write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);\n",
      " 261 :       vars.resize(vars_vec.size());\n",
      " 262 :       for (int i = 0; i < vars.size(); ++i)\n",
      " 263 :         vars(i) = vars_vec[i];\n",
      " 264 :     }\n",
      " 265 : \n",
      " 266 :     std::string model_name() const {\n",
      " 267 :         return \"modele188e6315_16a540c6086086816528e4524def24d9\";\n",
      " 268 :     }\n",
      " 269 : \n",
      " 270 : \n",
      " 271 :     void constrained_param_names(std::vector<std::string>& param_names__,\n",
      " 272 :                                  bool include_tparams__ = true,\n",
      " 273 :                                  bool include_gqs__ = true) const {\n",
      " 274 :         std::stringstream param_name_stream__;\n",
      " 275 :         param_name_stream__.str(std::string());\n",
      " 276 :         param_name_stream__ << \"y\";\n",
      " 277 :         param_names__.push_back(param_name_stream__.str());\n",
      " 278 : \n",
      " 279 :         if (!include_gqs__ && !include_tparams__) return;\n",
      " 280 : \n",
      " 281 :         if (include_tparams__) {\n",
      " 282 :         }\n",
      " 283 : \n",
      " 284 :         if (!include_gqs__) return;\n",
      " 285 :     }\n",
      " 286 : \n",
      " 287 : \n",
      " 288 :     void unconstrained_param_names(std::vector<std::string>& param_names__,\n",
      " 289 :                                    bool include_tparams__ = true,\n",
      " 290 :                                    bool include_gqs__ = true) const {\n",
      " 291 :         std::stringstream param_name_stream__;\n",
      " 292 :         param_name_stream__.str(std::string());\n",
      " 293 :         param_name_stream__ << \"y\";\n",
      " 294 :         param_names__.push_back(param_name_stream__.str());\n",
      " 295 : \n",
      " 296 :         if (!include_gqs__ && !include_tparams__) return;\n",
      " 297 : \n",
      " 298 :         if (include_tparams__) {\n",
      " 299 :         }\n",
      " 300 : \n",
      " 301 :         if (!include_gqs__) return;\n",
      " 302 :     }\n",
      " 303 : \n",
      " 304 : }; // model\n",
      " 305 : \n",
      " 306 : }  // namespace\n",
      " 307 : \n",
      " 308 : typedef modele188e6315_16a540c6086086816528e4524def24d9_namespace::modele188e6315_16a540c6086086816528e4524def24d9 stan_model;\n",
      " 309 : \n",
      " 310 : #ifndef USING_R\n",
      " 311 : \n",
      " 312 : stan::model::model_base& new_model(\n",
      " 313 :         stan::io::var_context& data_context,\n",
      " 314 :         unsigned int seed,\n",
      " 315 :         std::ostream* msg_stream) {\n",
      " 316 :   stan_model* m = new stan_model(data_context, seed, msg_stream);\n",
      " 317 :   return *m;\n",
      " 318 : }\n",
      " 319 : \n",
      " 320 : #endif\n",
      " 321 : \n",
      " 322 : \n",
      " 323 : \n",
      " 324 : #include <rstan_next/stan_fit.hpp>\n",
      " 325 : \n",
      " 326 : struct stan_model_holder {\n",
      " 327 :     stan_model_holder(rstan::io::rlist_ref_var_context rcontext,\n",
      " 328 :                       unsigned int random_seed)\n",
      " 329 :     : rcontext_(rcontext), random_seed_(random_seed)\n",
      " 330 :      {\n",
      " 331 :      }\n",
      " 332 : \n",
      " 333 :    //stan::math::ChainableStack ad_stack;\n",
      " 334 :    rstan::io::rlist_ref_var_context rcontext_;\n",
      " 335 :    unsigned int random_seed_;\n",
      " 336 : };\n",
      " 337 : \n",
      " 338 : Rcpp::XPtr<stan::model::model_base> model_ptr(stan_model_holder* smh) {\n",
      " 339 :   Rcpp::XPtr<stan::model::model_base> model_instance(new stan_model(smh->rcontext_, smh->random_seed_), true);\n",
      " 340 :   return model_instance;\n",
      " 341 : }\n",
      " 342 : \n",
      " 343 : Rcpp::XPtr<rstan::stan_fit_base> fit_ptr(stan_model_holder* smh) {\n",
      " 344 :   return Rcpp::XPtr<rstan::stan_fit_base>(new rstan::stan_fit(model_ptr(smh), smh->random_seed_), true);\n",
      " 345 : }\n",
      " 346 : \n",
      " 347 : std::string model_name(stan_model_holder* smh) {\n",
      " 348 :   return model_ptr(smh).get()->model_name();\n",
      " 349 : }\n",
      " 350 : \n",
      " 351 : RCPP_MODULE(stan_fit4modele188e6315_16a540c6086086816528e4524def24d9_mod){\n",
      " 352 :   Rcpp::class_<stan_model_holder>(\"stan_fit4modele188e6315_16a540c6086086816528e4524def24d9\")\n",
      " 353 :   .constructor<rstan::io::rlist_ref_var_context, unsigned int>()\n",
      " 354 :   .method(\"model_ptr\", &model_ptr)\n",
      " 355 :   .method(\"fit_ptr\", &fit_ptr)\n",
      " 356 :   .method(\"model_name\", &model_name)\n",
      " 357 :   ;\n",
      " 358 : }\n",
      " 359 : \n",
      " 360 : \n",
      " 361 : // declarations\n",
      " 362 : extern \"C\" {\n",
      " 363 : SEXP filee619bd1f9( ) ;\n",
      " 364 : }\n",
      " 365 : \n",
      " 366 : // definition\n",
      " 367 : SEXP filee619bd1f9() {\n",
      " 368 :  return Rcpp::wrap(\"16a540c6086086816528e4524def24d9\");\n",
      " 369 : }\n",
      "\n",
      "stn_md> fit <- sampling(mod, data = list(y_mean = 0))\n",
      "\n",
      "SAMPLING FOR MODEL '16a540c6086086816528e4524def24d9' NOW (CHAIN 1).\n",
      "Chain 1: \n",
      "Chain 1: Gradient evaluation took 1e-05 seconds\n",
      "Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.1 seconds.\n",
      "Chain 1: Adjust your expectations accordingly!\n",
      "Chain 1: \n",
      "Chain 1: \n",
      "Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)\n",
      "Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)\n",
      "Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)\n",
      "Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)\n",
      "Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)\n",
      "Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)\n",
      "Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)\n",
      "Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)\n",
      "Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)\n",
      "Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)\n",
      "Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)\n",
      "Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)\n",
      "Chain 1: \n",
      "Chain 1:  Elapsed Time: 0.008006 seconds (Warm-up)\n",
      "Chain 1:                0.00728 seconds (Sampling)\n",
      "Chain 1:                0.015286 seconds (Total)\n",
      "Chain 1: \n",
      "\n",
      "SAMPLING FOR MODEL '16a540c6086086816528e4524def24d9' NOW (CHAIN 2).\n",
      "Chain 2: \n",
      "Chain 2: Gradient evaluation took 6e-06 seconds\n",
      "Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.06 seconds.\n",
      "Chain 2: Adjust your expectations accordingly!\n",
      "Chain 2: \n",
      "Chain 2: \n",
      "Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)\n",
      "Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)\n",
      "Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)\n",
      "Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)\n",
      "Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)\n",
      "Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)\n",
      "Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)\n",
      "Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)\n",
      "Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)\n",
      "Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)\n",
      "Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)\n",
      "Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)\n",
      "Chain 2: \n",
      "Chain 2:  Elapsed Time: 0.007933 seconds (Warm-up)\n",
      "Chain 2:                0.008112 seconds (Sampling)\n",
      "Chain 2:                0.016045 seconds (Total)\n",
      "Chain 2: \n",
      "\n",
      "SAMPLING FOR MODEL '16a540c6086086816528e4524def24d9' NOW (CHAIN 3).\n",
      "Chain 3: \n",
      "Chain 3: Gradient evaluation took 6e-06 seconds\n",
      "Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.06 seconds.\n",
      "Chain 3: Adjust your expectations accordingly!\n",
      "Chain 3: \n",
      "Chain 3: \n",
      "Chain 3: Iteration:    1 / 2000 [  0%]  (Warmup)\n",
      "Chain 3: Iteration:  200 / 2000 [ 10%]  (Warmup)\n",
      "Chain 3: Iteration:  400 / 2000 [ 20%]  (Warmup)\n",
      "Chain 3: Iteration:  600 / 2000 [ 30%]  (Warmup)\n",
      "Chain 3: Iteration:  800 / 2000 [ 40%]  (Warmup)\n",
      "Chain 3: Iteration: 1000 / 2000 [ 50%]  (Warmup)\n",
      "Chain 3: Iteration: 1001 / 2000 [ 50%]  (Sampling)\n",
      "Chain 3: Iteration: 1200 / 2000 [ 60%]  (Sampling)\n",
      "Chain 3: Iteration: 1400 / 2000 [ 70%]  (Sampling)\n",
      "Chain 3: Iteration: 1600 / 2000 [ 80%]  (Sampling)\n",
      "Chain 3: Iteration: 1800 / 2000 [ 90%]  (Sampling)\n",
      "Chain 3: Iteration: 2000 / 2000 [100%]  (Sampling)\n",
      "Chain 3: \n",
      "Chain 3:  Elapsed Time: 0.007976 seconds (Warm-up)\n",
      "Chain 3:                0.007257 seconds (Sampling)\n",
      "Chain 3:                0.015233 seconds (Total)\n",
      "Chain 3: \n",
      "\n",
      "SAMPLING FOR MODEL '16a540c6086086816528e4524def24d9' NOW (CHAIN 4).\n",
      "Chain 4: \n",
      "Chain 4: Gradient evaluation took 6e-06 seconds\n",
      "Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.06 seconds.\n",
      "Chain 4: Adjust your expectations accordingly!\n",
      "Chain 4: \n",
      "Chain 4: \n",
      "Chain 4: Iteration:    1 / 2000 [  0%]  (Warmup)\n",
      "Chain 4: Iteration:  200 / 2000 [ 10%]  (Warmup)\n",
      "Chain 4: Iteration:  400 / 2000 [ 20%]  (Warmup)\n",
      "Chain 4: Iteration:  600 / 2000 [ 30%]  (Warmup)\n",
      "Chain 4: Iteration:  800 / 2000 [ 40%]  (Warmup)\n",
      "Chain 4: Iteration: 1000 / 2000 [ 50%]  (Warmup)\n",
      "Chain 4: Iteration: 1001 / 2000 [ 50%]  (Sampling)\n",
      "Chain 4: Iteration: 1200 / 2000 [ 60%]  (Sampling)\n",
      "Chain 4: Iteration: 1400 / 2000 [ 70%]  (Sampling)\n",
      "Chain 4: Iteration: 1600 / 2000 [ 80%]  (Sampling)\n",
      "Chain 4: Iteration: 1800 / 2000 [ 90%]  (Sampling)\n",
      "Chain 4: Iteration: 2000 / 2000 [100%]  (Sampling)\n",
      "Chain 4: \n",
      "Chain 4:  Elapsed Time: 0.00812 seconds (Warm-up)\n",
      "Chain 4:                0.007255 seconds (Sampling)\n",
      "Chain 4:                0.015375 seconds (Total)\n",
      "Chain 4: \n",
      "\n",
      "stn_md> fit2 <- sampling(mod, data = list(y_mean = 5))\n",
      "\n",
      "SAMPLING FOR MODEL '16a540c6086086816528e4524def24d9' NOW (CHAIN 1).\n",
      "Chain 1: \n",
      "Chain 1: Gradient evaluation took 8e-06 seconds\n",
      "Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.08 seconds.\n",
      "Chain 1: Adjust your expectations accordingly!\n",
      "Chain 1: \n",
      "Chain 1: \n",
      "Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)\n",
      "Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)\n",
      "Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)\n",
      "Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)\n",
      "Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)\n",
      "Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)\n",
      "Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)\n",
      "Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)\n",
      "Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)\n",
      "Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)\n",
      "Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)\n",
      "Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)\n",
      "Chain 1: \n",
      "Chain 1:  Elapsed Time: 0.009345 seconds (Warm-up)\n",
      "Chain 1:                0.008228 seconds (Sampling)\n",
      "Chain 1:                0.017573 seconds (Total)\n",
      "Chain 1: \n",
      "\n",
      "SAMPLING FOR MODEL '16a540c6086086816528e4524def24d9' NOW (CHAIN 2).\n",
      "Chain 2: \n",
      "Chain 2: Gradient evaluation took 6e-06 seconds\n",
      "Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.06 seconds.\n",
      "Chain 2: Adjust your expectations accordingly!\n",
      "Chain 2: \n",
      "Chain 2: \n",
      "Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)\n",
      "Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)\n",
      "Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)\n",
      "Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)\n",
      "Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)\n",
      "Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)\n",
      "Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)\n",
      "Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)\n",
      "Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)\n",
      "Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)\n",
      "Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)\n",
      "Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)\n",
      "Chain 2: \n",
      "Chain 2:  Elapsed Time: 0.008202 seconds (Warm-up)\n",
      "Chain 2:                0.008758 seconds (Sampling)\n",
      "Chain 2:                0.01696 seconds (Total)\n",
      "Chain 2: \n",
      "\n",
      "SAMPLING FOR MODEL '16a540c6086086816528e4524def24d9' NOW (CHAIN 3).\n",
      "Chain 3: \n",
      "Chain 3: Gradient evaluation took 5e-06 seconds\n",
      "Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.05 seconds.\n",
      "Chain 3: Adjust your expectations accordingly!\n",
      "Chain 3: \n",
      "Chain 3: \n",
      "Chain 3: Iteration:    1 / 2000 [  0%]  (Warmup)\n",
      "Chain 3: Iteration:  200 / 2000 [ 10%]  (Warmup)\n",
      "Chain 3: Iteration:  400 / 2000 [ 20%]  (Warmup)\n",
      "Chain 3: Iteration:  600 / 2000 [ 30%]  (Warmup)\n",
      "Chain 3: Iteration:  800 / 2000 [ 40%]  (Warmup)\n",
      "Chain 3: Iteration: 1000 / 2000 [ 50%]  (Warmup)\n",
      "Chain 3: Iteration: 1001 / 2000 [ 50%]  (Sampling)\n",
      "Chain 3: Iteration: 1200 / 2000 [ 60%]  (Sampling)\n",
      "Chain 3: Iteration: 1400 / 2000 [ 70%]  (Sampling)\n",
      "Chain 3: Iteration: 1600 / 2000 [ 80%]  (Sampling)\n",
      "Chain 3: Iteration: 1800 / 2000 [ 90%]  (Sampling)\n",
      "Chain 3: Iteration: 2000 / 2000 [100%]  (Sampling)\n",
      "Chain 3: \n",
      "Chain 3:  Elapsed Time: 0.008119 seconds (Warm-up)\n",
      "Chain 3:                0.008086 seconds (Sampling)\n",
      "Chain 3:                0.016205 seconds (Total)\n",
      "Chain 3: \n",
      "\n",
      "SAMPLING FOR MODEL '16a540c6086086816528e4524def24d9' NOW (CHAIN 4).\n",
      "Chain 4: \n",
      "Chain 4: Gradient evaluation took 5e-06 seconds\n",
      "Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.05 seconds.\n",
      "Chain 4: Adjust your expectations accordingly!\n",
      "Chain 4: \n",
      "Chain 4: \n",
      "Chain 4: Iteration:    1 / 2000 [  0%]  (Warmup)\n",
      "Chain 4: Iteration:  200 / 2000 [ 10%]  (Warmup)\n",
      "Chain 4: Iteration:  400 / 2000 [ 20%]  (Warmup)\n",
      "Chain 4: Iteration:  600 / 2000 [ 30%]  (Warmup)\n",
      "Chain 4: Iteration:  800 / 2000 [ 40%]  (Warmup)\n",
      "Chain 4: Iteration: 1000 / 2000 [ 50%]  (Warmup)\n",
      "Chain 4: Iteration: 1001 / 2000 [ 50%]  (Sampling)\n",
      "Chain 4: Iteration: 1200 / 2000 [ 60%]  (Sampling)\n",
      "Chain 4: Iteration: 1400 / 2000 [ 70%]  (Sampling)\n",
      "Chain 4: Iteration: 1600 / 2000 [ 80%]  (Sampling)\n",
      "Chain 4: Iteration: 1800 / 2000 [ 90%]  (Sampling)\n",
      "Chain 4: Iteration: 2000 / 2000 [100%]  (Sampling)\n",
      "Chain 4: \n",
      "Chain 4:  Elapsed Time: 0.007965 seconds (Warm-up)\n",
      "Chain 4:                0.007336 seconds (Sampling)\n",
      "Chain 4:                0.015301 seconds (Total)\n",
      "Chain 4: \n"
     ]
    }
   ],
   "source": [
    "example(stan_model, package = \"rstan\", run.dontrun = TRUE)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "id": "1fdece74",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-07-30T03:01:17.970191Z",
     "iopub.status.busy": "2022-07-30T03:01:17.967544Z",
     "iopub.status.idle": "2022-07-30T03:01:18.000164Z",
     "shell.execute_reply": "2022-07-30T03:01:17.998170Z"
    },
    "papermill": {
     "duration": 0.052695,
     "end_time": "2022-07-30T03:01:18.003156",
     "exception": false,
     "start_time": "2022-07-30T03:01:17.950461",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [
    {
     "data": {
      "text/html": [
       "4"
      ],
      "text/latex": [
       "4"
      ],
      "text/markdown": [
       "4"
      ],
      "text/plain": [
       "[1] 4"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "parallel::detectCores()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "id": "16b9d1f0",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-07-30T03:01:18.014517Z",
     "iopub.status.busy": "2022-07-30T03:01:18.012859Z",
     "iopub.status.idle": "2022-07-30T03:02:16.337105Z",
     "shell.execute_reply": "2022-07-30T03:02:16.334998Z"
    },
    "papermill": {
     "duration": 58.333068,
     "end_time": "2022-07-30T03:02:16.339898",
     "exception": false,
     "start_time": "2022-07-30T03:01:18.006830",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Installing package into ‘/usr/local/lib/R/site-library’\n",
      "(as ‘lib’ is unspecified)\n",
      "\n"
     ]
    }
   ],
   "source": [
    "install.packages(\"brms\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "id": "c6cb4976",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-07-30T03:02:16.351180Z",
     "iopub.status.busy": "2022-07-30T03:02:16.349415Z",
     "iopub.status.idle": "2022-07-30T03:02:18.195247Z",
     "shell.execute_reply": "2022-07-30T03:02:18.192955Z"
    },
    "papermill": {
     "duration": 1.854139,
     "end_time": "2022-07-30T03:02:18.197925",
     "exception": false,
     "start_time": "2022-07-30T03:02:16.343786",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Loading required package: Rcpp\n",
      "\n",
      "Loading 'brms' package (version 2.17.0). Useful instructions\n",
      "can be found by typing help('brms'). A more detailed introduction\n",
      "to the package is available through vignette('brms_overview').\n",
      "\n",
      "\n",
      "Attaching package: ‘brms’\n",
      "\n",
      "\n",
      "The following object is masked from ‘package:rstan’:\n",
      "\n",
      "    loo\n",
      "\n",
      "\n",
      "The following object is masked from ‘package:stats’:\n",
      "\n",
      "    ar\n",
      "\n",
      "\n"
     ]
    }
   ],
   "source": [
    "library(brms)"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "R",
   "language": "R",
   "name": "ir"
  },
  "language_info": {
   "codemirror_mode": "r",
   "file_extension": ".r",
   "mimetype": "text/x-r-source",
   "name": "R",
   "pygments_lexer": "r",
   "version": "4.0.5"
  },
  "papermill": {
   "default_parameters": {},
   "duration": 965.419822,
   "end_time": "2022-07-30T03:02:18.323613",
   "environment_variables": {},
   "exception": null,
   "input_path": "__notebook__.ipynb",
   "output_path": "__notebook__.ipynb",
   "parameters": {},
   "start_time": "2022-07-30T02:46:12.903791",
   "version": "2.3.4"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
