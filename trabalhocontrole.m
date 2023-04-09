format long % aumento do numero de casas decimais
clear all;  % limpa o ambiente MatLab de todas as variaveis em uso
close all;  % fecha todas as janelas de plotagem que foram abertas
clc;        % Limpa a janela de comando do Matlab
s = tf('s');% Especificando a funcao de transferencia = s, que é a variavel de Laplace

%% Trabalho Computacional - Teoria de Controle I CEL038 Turma B
% Grupo: Bruno de Oliveira Melão
%        Josiano Aparecido Silva
%        Lucas Almeida Silva
%        Lucas Campos Mahia


% Item A: Root locus com reta de amortecimento 0.4
figure(1);
hold on;
% Gráfico do lugar das raízes
num=[7]; %Numerador
den=conv([1 1.5 0],[1 2.5]); %Denominador
rlocus(num, den); 
v=[-5 5 -5 5];  % Escala do gráfico
axis(v); % Formato do gráfico
sgrid;
title('Gráfico do lugar das raízes do sistema não compensado'); %Título 
xlabel('Eixo Real'); %Legenda do eixo x
ylabel('Eixo Imaginário'); %Legenda do eixo y
 
%Determinação da Função de Transferência de malha fechada e demarcação dos
%polos de malha fechada no Root Locus (*)
FTMA= 7/(s*(s+1.5)*(s+2.5)); 
FTMF= feedback(FTMA,1); %Função de transferência de malha fechada através da função feedback
[p,~]=rlocus(FTMA, 1);% Armazenamento dos polos de malha fechada para ganho unitário
wn1=abs(p(2));
zeta1=abs(real(p(2)))/wn1;
sgrid([0.4 zeta1],[2 wn1]); %Reta de amortecimento para zeta=0.4 e círculo de wn de 2 rad/s e para os polos dominantes
plot(p,'*'); % Polos de malha fechada para ganho=1 marcados no gráfico com *

%Item B: Gráfico da resposta no tempo de c(t) quando aplicado um degrau
%unitário na entrada

figure(2);
hold on;

t=0:0.01:50; % Vetor de tempo de 0 a 50 segundos com intervalos de 0.01
c1=step(FTMF,t); %Função step para resposta ao degrau
plot(t,c1,'-');
grid;
title('Resposta ao degrau unitário do sistema não compensado'); %Título do gráfico da fig2
xlabel('Tempo(segundos)'); % Legenda eixo x da fig2
ylabel('Amplitude'); % Legenda eixo y da fig2

%Item C: Determinar máximo sobressinal percentual, f, wd. Comparar wd com o
%obtido pelos dados do item A

[m,t1]=max(c1); % m=máximo sobressinal e tp=tempo onde ocorreu o máximo
tempo_pico=(t1-1)*0.01; %Fórmula para o tempo de pico 
[minimo,t2]=min(c1(t1:end)); %Determinação do mínimo depois do tempo de pico 
tempo_minimo=(t2+t1)*0.01;
periodo=2*(tempo_minimo-tempo_pico); % Do tempo de pico ao próximo tempo de mínimo temos meio periodo
freq=1/periodo; %Fórmula para cálculo da frequência
wd2=2*pi*freq; %Fórmula para frequência amortecida
Mp=(m-1)*100; % Máximo sobressinal percentual 

%Comparação entre wd1 e wd2
% wd1 é a parte imaginária do p(2)= 1.391941090707505
% wd2= 1.383961521405195
% Como esperado os valores de wd1 e wd2 são aproximadamente iguais, pois o
% gráfico da resposta ao degrau deve ter uma frequência wd igual a parte
% imaginária do polo complexo de malha fechada para ganho unitário no sistema
% de terceira ordem.


%Item D
figure(3);
hold on;

t=0:0.01:200;% Vetor de tempo com 200seg de simulação com intervalos de 0.01segundos
c2=lsim(FTMF,t,t); %Função lsim para resposta a rampa
plot(t,c2,'-');
grid;
title('Resposta a rampa no sistema não compensado'); %Titulo do gráfico da figura 3
xlabel('Tempo(segundos)'); % Legenda eixo x da fig3
ylabel('Amplitude'); % Legenda eixo y da fig3
plot(t,t); % Rampa unitária como referência

%Valor da constante de velocidade Kv e erro de regime permanente
syms S %criando váriavel simbolica para calcular limite
f=7/((S+1.5)*(S+2.5)); %Função de transferência utilizando variável S
Kv=double(limit(f,S,0)); %Cálculo do limite quando S tende a 0 determinando o erro de velocidade
ess=1/Kv; %Erro de regime permanente
% Valores obtidos
%Kv= 1.866666666666667
%ess= 0.535714285714286

erro_rampa=(max(t)-max(c2)); % Erro em regime permanente apresentado pelo gráfico de saída
%erro_rampa= 0.535714285713709

% O erro obtido pelo gráfico da resposta a rampa é igual ao calculado,
% confirmando que essa é a resposta a rampa do sistema.


%Item E: Projetar um compensador 
% Parâmetros a serem considerados:
%   zeta=0.4;
%   wn=2;
%   ess=ess'/5
%   Mp<0.36
%   ts<=3.8s

%Encontrar polo desejado a partir do zeta=0.4 e wn=2
% p=-sigma +- j*wd_itemA
zeta=0.4;
wn=2;
sigma=zeta*wn;
wd=wn*sqrt(1-(zeta^2));
polosDes=[-sigma+1i*wd -sigma-1i*wd];

%Para o cálculo da deficiência angular, é preciso saber quais os polos de
%malha aberta:
polosMA=[0 -1.5 -2.5]; %polos de malha aberta

%Cálculo da deficiência angular
phi=pi-(angle(polosDes(1)-polosMA(1)))-(angle(polosDes(1)-polosMA(2)))-(angle(polosDes(1)-polosMA(3)));
Cont_angular=-phi; % O compensador deve contribuir com -phi para que o sistema seja compensado

%Parte para cálculo do polo do compensador por avanço de fase 

aux=1; % Variavel auxiliar para vetor de resultados
cont=1; % Variavel auxiliar para vetor do compensador
resultados=[]; % Declaração do vetor de resultados para receber os valores relacionados ao atraso 
compensador=[]; % Declaração do vetor de resultados para receber os valores relacionados ao compensador como um todo

for zav = -0.5: -0.1: -2 %Zero de avanço escolhido na faixa de -0.5 a -2 com intervalos de 0.1
 
    angpz=angle(polosDes(1)-zav);% Cálculo do ângulo entre o polo e zero
    beta_ang=angpz-Cont_angular; %beta_ang: ângulo entre o polo de avanço e polo desejado
    %Pela relação entre os ângulos do triângulo temos que 
    % Pavanço= imag(polo desejado)/tan(angulo entre pAvanço e polo
    % desejado)+abs(real(polo desejado))
    pav=-(wd/(abs(tan(beta_ang)))+sigma); % Como utilizei zav negativo, pav também é negativo
    
    % Pela condição de módulo temos que 
    % abs(Kc*Gavanco*G(s))=1
    % Determinamos então o Kc
    Kc = 1/(abs((polosDes(1)-zav)/(polosDes(1)-pav))*abs(7/(polosDes(1)*(polosDes(1)+1.5)*(polosDes(1)+2.5))));
    Gav=Kc*(s-zav)/(s-pav); %Função de transferência do compensador avanço
    
    % Como s+1/T1=0 quando s=zav, isolando T1 temos:
    T1=1/(-zav); 
    % Como s+gama/T1=0 quando s=pav, isolando gama temos:
    gama=-pav*T1;
    
    %Projeto do compensador atraso
    %Kn=Kv/5
    %Kn=Kv*Kc/gama*beta
    beta=5*gama/Kc; %Fórmula do beta
    
    %Escolher T2 de modo que a condição de ângulo e modulo sejam
    %satisfeitas
    
    for T2=1:1:10 %Varredura do T2
        anguloAtraso=-rad2deg(angle(((polosDes(1))+(1/T2))/((polosDes(1))+(1/(beta*T2))))); %Cálculo do ângulo do compensador atraso
        moduloAtraso=abs((((polosDes(1))+(1/T2))/((polosDes(1))+(1/(beta*T2))))); %Cálculo do módulo do compensador atraso
     
        if anguloAtraso>0 && anguloAtraso<=5 && moduloAtraso>=0.96 %Condições de módulo e ângulo
            Gatraso=(s+1/T2)/(s+1/(beta*T2)); %Função de transferência do compensador atraso
            Gc=Gav*Gatraso; %Função de transferência do compensador avanço-atraso
            
            FTMAcompensado=FTMA*Gc; %Função de transferência de malha aberta do sistema compensado
            FTMFcompensado=feedback(FTMAcompensado,1); %Função de transferência de malha aberta do sistema compensado
            
            t=0:0.01:50;% Tempo da simulação para entrada em degrau
            rdegrau=step(FTMFcompensado,t);%Com a FTMF compensada, aplica-se um degrau na entrada pela função step
            
            m2=max(rdegrau); %Posição do máximo sobressinal da resposta ao degrau
            MP2=(m2-1)*100; %Máximo sobressinal percentual
            
            if m2<1.36 %Condição do projeto para sobressinal de até 36% acima
                
                %Busca de trás para frente para fazer menor esforço
                %computacional
                %Para isso, definimos o comprimento do vetor de tempo e
                %fazemos a busca de trás para frente
                comp=length(t); %função length
                
                % Busca inversa pela faixa 5%
                while rdegrau(comp)>0.95 &&rdegrau(comp)<1.05
                    comp=comp-1;
                end
                %Através da posição de comp se calcula o tempo de
                %assentamento
                tassentamento=(comp-1)*0.01 ; %Posição no vetor multiplicado pelo tamanho dos intervalos
                if tassentamento <= 3.8 %Condição do tempo de estabilização
                    resultados(aux,:)= [T2 m2 tassentamento]; %Armazenamento na matriz resultados dos valores de T2, m2 e tassentamento que satisfazem as restrições
                    
                    aux=aux+1; %iteração da variável auxiliar
                end
            end
        end
    end
     if length(resultados)>0 % Verificar se houve algum resultado que satisfez as condições (importante porque ainda está dentro de um loop)
        resultados_classificados = sortrows(resultados, 2); %Classificar os sobressinais(coluna 2) de modo crescente
        %Armazenamento na matriz compensador os valores relacionados ao
        %compensador que são opções
        compensador(cont,:) = [angpz pav zav Kc T1 gama beta resultados_classificados(1,1) resultados_classificados(1,2) resultados_classificados(1,3)]; 
                      % Ordem [angpz pav zav Kc T1 gama beta T2 menorsobressinal tassentamento]
        cont=cont+ 1;
     end

end
%Agora definimos o melhor compensador avanço e atraso depois de ter
%determinado o melhor atraso para cada avanço 
%Como a matriz está ordenada por menor sobressinal, foi escolhido a linha
% 1 da matriz como escolha de projeto por ter o melhor sobressinal

%Escolha dos parâmetros do compensador
Kc_escolhido = compensador(1, 4);
T1_escolhido = compensador(1, 5);
T2_escolhido = compensador(1, 8);
gama_escolhido = compensador(1, 6);
beta_escolhido = compensador(1, 7);

Gav_escolhido=Kc_escolhido*((s+(1/T1_escolhido))/(s+(gama_escolhido/T1_escolhido))); %Função de transferência do compensador avanço escolhido
Gatraso_escolhido=((s+(1/T2_escolhido))/(s+(1/(beta_escolhido*T2_escolhido)))); %Função de transferência do compensador atraso escolhido
Gc_escolhido=Gav_escolhido*Gatraso_escolhido;%Função de transferência do compensador avanço-atraso escolhido

FTMA_escolhido=FTMA*Gc_escolhido; %Função de transferência de malha aberta do sistema compensado escolhido
FTMF_escolhido=feedback(FTMA_escolhido,1);%Função de transferência de malha fechada do compensador escolhido

%Root Locus do sistema compensado com polos dominantes de ganho unitário
%marcado com * e zeta e wn especificados
figure(4);
hold on;
rlocus(FTMA_escolhido);
[pcomp,~]=rlocus(FTMA_escolhido,1);
title('Grafico do Lugar das Raizes - Sistema Compensado');
xlabel('Eixo Real');
ylabel('Eixo Imaginario');
sgrid([0.4],[2]);
plot(pcomp,'*');

%Resposta ao degrau do sistema compensado com sobressinal e tempo de
%assentamento de acordo com projeto
figure(5)
hold on;
t=0:0.01:50;
c_compensado=step(FTMF_escolhido,t); 
plot(t,c_compensado,'-');
grid;
% Definição do título do gráfico e do eixos
title('Grafico da resposta ao degrau dos sistema compensado');
xlabel('Tempo');
ylabel('Amplitude');

% Resposta a rampa unitária do sistema compensado com erro de regime
% permanente 5 vezes menor que o erro calculado e obtido no gráfico do item
% D
figure(6)
hold on;
t=0:0.1:200;
c_rampa=lsim(FTMF_escolhido,t,t); % Aplicação da rampa à FT de malha fechada do sistema não comoensado no intervalo de tempo t2
plot(t,c_rampa,'-');
plot(t,t);
grid;
% Definição do título do gráfico e do eixos
title('Grafico da resposta a rampa unitaria dos sistema compensado');
xlabel('Tempo');
ylabel('Amplitude');
erro_rampaComp=(max(t)-max(c_rampa)); %Verificando que erro da rampa diminuiu 5 vezes
%Novo erro da entrada a rampa em regime permanente
%ess=0.107142857525417;