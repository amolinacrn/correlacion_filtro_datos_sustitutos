\chapter{Metodologia}
%não se esqueça de definir uma label única para utilizar no comando \ref
\label{chap:metodologia}
Nesta pesquisa, analisamos a atividade neuronal a partir de dados eletrofisiológicos. O conjunto de dados consiste de três horas de registros eletrofisiológicos. Para o cálculo das correlações, usamos a função de correlação (Equação. \ref{ecuacion:NCCH}) descrita no capítulo \ref{Cap:correlations}. Para este estudo, dividimos os bancos de dados em segmentos de 250 segundos para minimizar ao máximo os efeitos negativos causados pela não estacionariedade das séries temporais. Nesses segmentos, os coeficientes de variação são calculados a cada 10 segundos e, em seguida, calculamos uma média para obter um único valor médio correspondente. Em seguida, para cada um desses intervalos de tempo, calculamos todas as correlações cruzadas entre pares de neurônios e construímos as matrizes de conectividade usando as técnicas descritas no capítulo \ref{Cap:correlations}. Depois, aplicamos um limiar para eliminar conexões ruidosas ou espúrias e, por fim, calculamos as métricas de rede mencionadas no capítulo \ref{redes_grafos}.

\section{Cálculo de correlações}

Uma das tarefas mais exaustivas neste trabalho é o cálculo de todas as correlações cruzadas, esta etapa é uma das mais complicadas devido ao grande número de correlações que precisam ser calculadas. 
Mas não só isso, também é importante criar simultaneamente dados surrogados para poder aplicar o limiar que elimina correlações espúrias, portanto, o código deve calcular todas as correlações, gerar simultaneamente dados surrogados, calcular um limiar e avaliar as correlações se elas são significativas ou não. Aparentemente, todo este processo parece ser uma tarefa simples, porém, é um processo tedioso e computacionalmente intensivo. Para resolver este problema e ser capaz de calcular milhares de correlações, nós nos apoiamos no software Spicodyn \cite{spicodyn}, que é precisamente projetado para ter uma alta eficiência computacional no cálculo de correlações e na geração de dados surrogados. Spicodyn é um software que faz parte do Neuroimaging Tools and Resources Collaboratory (NITRC), que se dedica ao desenvolvimento de ambientes de conhecimento de fácil utilização, com foco na melhoria e adoção de ferramentas e recursos em neuroimagem, imagem ótica, neuroinformática clínica e neurociência computacional \cite{nitrc}.  Entretanto, há uma dificuldade que não é desprezível, Spicodyn é um software com um ambiente gráfico que é projetado e limitado a experimentos \textit{in vitro} \cite{invitro}, esta limitação impede que um pesquisador possa usar o sistema para analisar conjuntos de dados pertencentes a um número arbitrário de neurônios que foram detectados em experimentos \textit{in vivo} \cite{invivo}. \\

\begin{figure}[ht!]
\centering
\caption{\textmd{Variáveis locais a serem definidas antes de iniciar o processo de cálculo da correlação. As linhas 8 a 13 definem a janela de correlação (\textit{window}), o número de dados surrogados (\textit{num-surr}), a freqüência de amostragem dos dados (\textit{fs}), a taxa mínima de disparo (\textit{mfr}), o tamanho do bin (\textit{bin}) e a janela de shufling (\textit{wsurrgte}), a variável \textit{set-exp} define todos os experimentos que são analisados \{ExpMar07, ExpJan29, ExpDez20, ExpMar04, ExpJan14, ExpMar10, ExpJan21\}}}
\label{fig:variable}
\fcolorbox{gray}{white}{\includegraphics[width=0.97\textwidth]{images/captitulo1/variables.png}
%\par\medskip\ABNTEXfontereduzida\selectfont\textbf{Fonte:}{
}
 \par\medskip\ABNTEXfontereduzida\selectfont\textbf{Fonte:} Elaborada pelo autor (2023) \par\medskip
\end{figure}


Em experimentação \textit{in vitro}, os eletrodos são colocados em uma geometria completamente diferente da forma como são distribuídos em uma sonda de silício. Além disso, os dados aqui analisados são classificados por disparos (Spicodyn não é projetado para fazer a classificação por disparos), o que significa que o número de trens de disparos presentes é completamente arbitrário e NÃO está limitado ao número de eletrodos. A classificação de disparos é muito importante porque é uma garantia de que os dados obtidos após a classificação podem ser identificados para cada um dos diferentes neurônios individuais (dados SUA). Para fazer uso desta importante ferramenta computacional, foi necessário baixar o sistema com todos os seus pacotes, extrair os algoritmos necessários e depois juntá-los um a um, ou seja, o sistema é composto de dezenas de funções escritas em milhares de linhas de código, portanto, era inevitável ter que chegar aos scripts fonte, identificar os métodos e atributos que estamos procurando, extraí-los e finalmente reconstruí-los para que sejam totalmente funcionais a qualquer conjunto de dados. Este sistema é escrito em uma linguagem de programação C\# e o pacote  resultante contém aproximadamente 2.000 linhas de código.\\

Para a construção das matrizes de conectividade funcional, seguimos um processo que envolve várias etapas. Inicialmente, o coeficiente de variação ($CV$) é calculado em intervalos de tempo de 10 segundos. Esses resultados são organizados e, posteriormente, agrupados em conjuntos de CVs adjacentes que são calculados como média. Portanto, para cada valor $\langle CV \rangle$ específico, os dados correspondentes são extraídos. Isso resulta em séries temporais que abrangem um período de 250 segundos com um nível específico de variabilidade da atividade neural. Esse processo produz trens de disparos caracterizados por variabilidade neuronal específica, permitindo a geração subsequente de matrizes de conectividade funcional.\\

\begin{figure}[ht!]
\centering

\caption{\textmd{Conjunto de correlações cruzadas com sua respectiva matriz de conectividade com limiar. Aqui são mostrados apenas as correlações significativas. Os parâmetros usados para o cálculo são mostrados na figura \ref{fig:variable}. Os dados são do tipo SUA e correspondem ao experimento realizado no dia 7 de março do ano 2020. O gráfico é feito em python e o pacote \href{https://matplotlib.org/}{\texttt{Matplotlib}} 
%(https://doi.org/10.1103/PhysRevE.103.012415)
}}
\label{fig:_mcc_corr_}
\fcolorbox{gray}{white}{\includegraphics[width=1\textwidth]{images/captitulo1/mcc_corr.png}}
 \par\medskip\ABNTEXfontereduzida\selectfont\textbf{Fonte:} Elaborada pelo autor (2023)
%\par\medskip\ABNTEXfontereduzida\selectfont\textbf{Fonte:} \citeauthor{} (\citeyear{manualufpe2020}) \par\medskip
\end{figure}

Agora, para cada trem de disparos de 250 segundos, os dados são extraídos por neurônio. Esse agrupamento de disparos é importante porque permite o cálculo de todas as correlações de pares entre os diferentes neurônios. Após esse estágio, a próxima etapa é calcular as correlações com o pacote de código C\# que foi preparado para essa tarefa. O método que carrega todo o sistema tem todas as variáveis que precisam ser definidas. Na linha 4 da figura \ref{fig:variable} se definem os diretórios onde estão localizados os experimentos que se deseja analisar, aqui se analisam 6 experimentos diferentes \{ExpMar07, ExpJan29, ExpMar04, ExpJan14, ExpMar10, ExpJan21\}. As linhas 8 a 13 definem a janela de correlação (\textit{window}), o número de dados surrogados (\textit{num-surr}), a freqüência de amostragem dos dados (\textit{fs}), a taxa mínima de disparo (\textit{mfr}), o tamanho do bin (\textit{bin}) e a janela de shufling (\textit{wsurrgte}). \\

Observo que, durante o processo de divisão de dados, é possível que ocorram intervalos com um pequeno número de disparos. Nessa situação, é essencial que o código determine se a série temporal tem um número significativo de disparos ou não. Portanto, as correlações entre as várias séries temporais se tornarão efetivas se a expressão a seguir for satisfeita:
\begin{equation}
    f_{r}\leq\frac{N_{\textup{spike}}}{t_{\textup{spike}}/f_s}
\end{equation}

onde $N_{\textup{spike}}$, $t_{\textup{spike}}$, $f_s$, é o número total de disparos, o maior tempo de pico e a taxa de amostragem, da série temporal, respectivamente. $fr=0.1$, é a taxa mínima de disparo da resolução de código.\\

\begin{figure}[ht!]
\centering

\caption{\textmd{Coeficientes de variação para todos os experimentos. Cada $CV$ é calculado em segmentos ou janelas de tempo de 10 segundos , o tamanho do bin é de 50 milissegundos. O gráfico é feito em python e na biblioteca \href{https://matplotlib.org/}{\texttt{Matplotlib}} 
%(https://doi.org/10.1103/PhysRevE.103.012415)
}}
\label{fig:valores_de_cv}
\fcolorbox{gray}{white}{\includegraphics[width=1\textwidth]{images/captitulo1/cv.png}}
 \par\medskip\ABNTEXfontereduzida\selectfont\textbf{Fonte:} Elaborada pelo autor (2023)
%\par\medskip\ABNTEXfontereduzida\selectfont\textbf{Fonte:} \citeauthor{} (\citeyear{manualufpe2020}) \par\medskip
\end{figure}

O resultado depois de todo o processo computacional,  é o conjunto de correlações cruzadas que precisam ser organizados para construir a matriz de conectividade funcional (Fig. \ref{fig:_mcc_corr_}). Aqui a matriz de conectividade coincide com uma matriz de adjacência ponderada, mas isto nem sempre é o caso, em outras abordagens os resultados são uma matriz completa de tamanho $n \times n$, nestes casos, os valores na diagonal principal devem ser eliminados para depois obter um grafo não direcionado e sem loops \cite{fornito,Lee,Zampieri}.\\

Por outro lado, a Fig. \ref{fig:valores_de_cv} apresenta os valores de $CV$ calculados para cada experimento, enquanto a Fig. \ref{fig:cv_funcioal} apresenta cinco matrizes de conectividade funcional construídas para diferentes níveis médios de variabilidade que são destacados em cores (azul, verde, amarelo, magenta e vermelho), Observe como as entradas das matrizes de conectividade funcional se tornam mais fortes para cada valor médio de CV. Deve-se esclarecer que os pontos coloridos são aqueles cuja média é calculada para obter um único valor de $\langle CV \rangle$.

\begin{figure}[ht!]
\centering

\caption{\textmd{Coeficientes de variação para todos os experimentos. Cada $CV$ é calculado em segmentos ou janelas de tempo de 10 segundos , o tamanho do bin é de 50 milissegundos.. O gráfico é feito em python e na biblioteca \href{https://matplotlib.org/}{\texttt{Matplotlib}} 
%(https://doi.org/10.1103/PhysRevE.103.012415)
}}
\label{fig:cv_funcioal}
\fcolorbox{gray}{white}{\includegraphics[width=1\textwidth]{images/captitulo1/cv_funcioal.png}}
 \par\medskip\ABNTEXfontereduzida\selectfont\textbf{Fonte:} Elaborada pelo autor (2023)
%\par\medskip\ABNTEXfontereduzida\selectfont\textbf{Fonte:} \citeauthor{} (\citeyear{manualufpe2020}) \par\medskip
\end{figure}

Finalmente, quando todas as matrizes de conectividade conectadas estiverem prontas. Prosseguimos com o cálculo de todas as métricas apresentadas no capítulo 
 capítulo \ref{redes_grafos}. Para isso, são implementados códigos em Python. Nós nos focalizamos no cálculo do grau médio (Equação \ref{grado_medio_pesado}), densidade (Equação \ref{ecu:desidad}), Total de nós, nós desconectados e conectados, propensão ao mundo pequeno \ref{ecu:phi_swp}, caminhos mais curtos \ref{distancia_pesada}, agrupamento \ref{agrupacion_pesado} e modularidade. Mais informações sobre essa etapa podem ser encontradas na seção de resultados. 



