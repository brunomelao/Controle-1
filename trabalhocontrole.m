format long % aumento do numero de casas decimais
clear all;  % limpa o ambiente MatLab de todas as variaveis em uso
close all;  % fecha todas as janelas de plotagem que foram abertas
clc;        % Limpa a janela de comando do Matlab
s = tf('s');% Especificando a funcao de transferencia = s, que � a variavel de Laplace

%% Trabalho Computacional - Teoria de Controle I CEL038 Turma B
% Grupo: Bruno de Oliveira Mel�o
%        Josiano Aparecido Silva
%        Lucas Almeida Silva
%        Lucas Campos Mahia


% Item A: Root locus com reta de amortecimento 0.4
figure(1);
hold on;
% Gr�fico do lugar das ra�zes
num=[7]; %Numerador
den=conv([1 1.5 0],[1 2.5]); %Denominador
rlocus(num, den); 
v=[-5 5 -5 5];  % Escala do gr�fico
axis(v); % Formato do gr�fico
sgrid;
title('Gr�fico do lugar das ra�zes do sistema n�o compensado'); %T�tulo 
xlabel('Eixo Real'); %Legenda do eixo x
ylabel('Eixo Imagin�rio'); %Legenda do eixo y
 
%Determina��o da Fun��o de Transfer�ncia de malha fechada e demarca��o dos
%polos de malha fechada no Root Locus (*)
FTMA= 7/(s*(s+1.5)*(s+2.5)); 
FTMF= feedback(FTMA,1); %Fun��o de transfer�ncia de malha fechada atrav�s da fun��o feedback
[p,~]=rlocus(FTMA, 1);% Armazenamento dos polos de malha fechada para ganho unit�rio
wn1=abs(p(2));
zeta1=abs(real(p(2)))/wn1;
sgrid([0.4 zeta1],[2 wn1]); %Reta de amortecimento para zeta=0.4 e c�rculo de wn de 2 rad/s e para os polos dominantes
plot(p,'*'); % Polos de malha fechada para ganho=1 marcados no gr�fico com *

%Item B: Gr�fico da resposta no tempo de c(t) quando aplicado um degrau
%unit�rio na entrada

figure(2);
hold on;

t=0:0.01:50; % Vetor de tempo de 0 a 50 segundos com intervalos de 0.01
c1=step(FTMF,t); %Fun��o step para resposta ao degrau
plot(t,c1,'-');
grid;
title('Resposta ao degrau unit�rio do sistema n�o compensado'); %T�tulo do gr�fico da fig2
xlabel('Tempo(segundos)'); % Legenda eixo x da fig2
ylabel('Amplitude'); % Legenda eixo y da fig2

%Item C: Determinar m�ximo sobressinal percentual, f, wd. Comparar wd com o
%obtido pelos dados do item A

[m,t1]=max(c1); % m=m�ximo sobressinal e tp=tempo onde ocorreu o m�ximo
tempo_pico=(t1-1)*0.01; %F�rmula para o tempo de pico 
[minimo,t2]=min(c1(t1:end)); %Determina��o do m�nimo depois do tempo de pico 
tempo_minimo=(t2+t1)*0.01;
periodo=2*(tempo_minimo-tempo_pico); % Do tempo de pico ao pr�ximo tempo de m�nimo temos meio periodo
freq=1/periodo; %F�rmula para c�lculo da frequ�ncia
wd2=2*pi*freq; %F�rmula para frequ�ncia amortecida
Mp=(m-1)*100; % M�ximo sobressinal percentual 

%Compara��o entre wd1 e wd2
% wd1 � a parte imagin�ria do p(2)= 1.391941090707505
% wd2= 1.383961521405195
% Como esperado os valores de wd1 e wd2 s�o aproximadamente iguais, pois o
% gr�fico da resposta ao degrau deve ter uma frequ�ncia wd igual a parte
% imagin�ria do polo complexo de malha fechada para ganho unit�rio no sistema
% de terceira ordem.


%Item D
figure(3);
hold on;

t=0:0.01:200;% Vetor de tempo com 200seg de simula��o com intervalos de 0.01segundos
c2=lsim(FTMF,t,t); %Fun��o lsim para resposta a rampa
plot(t,c2,'-');
grid;
title('Resposta a rampa no sistema n�o compensado'); %Titulo do gr�fico da figura 3
xlabel('Tempo(segundos)'); % Legenda eixo x da fig3
ylabel('Amplitude'); % Legenda eixo y da fig3
plot(t,t); % Rampa unit�ria como refer�ncia

%Valor da constante de velocidade Kv e erro de regime permanente
syms S %criando v�riavel simbolica para calcular limite
f=7/((S+1.5)*(S+2.5)); %Fun��o de transfer�ncia utilizando vari�vel S
Kv=double(limit(f,S,0)); %C�lculo do limite quando S tende a 0 determinando o erro de velocidade
ess=1/Kv; %Erro de regime permanente
% Valores obtidos
%Kv= 1.866666666666667
%ess= 0.535714285714286

erro_rampa=(max(t)-max(c2)); % Erro em regime permanente apresentado pelo gr�fico de sa�da
%erro_rampa= 0.535714285713709

% O erro obtido pelo gr�fico da resposta a rampa � igual ao calculado,
% confirmando que essa � a resposta a rampa do sistema.


%Item E: Projetar um compensador 
% Par�metros a serem considerados:
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

%Para o c�lculo da defici�ncia angular, � preciso saber quais os polos de
%malha aberta:
polosMA=[0 -1.5 -2.5]; %polos de malha aberta

%C�lculo da defici�ncia angular
phi=pi-(angle(polosDes(1)-polosMA(1)))-(angle(polosDes(1)-polosMA(2)))-(angle(polosDes(1)-polosMA(3)));
Cont_angular=-phi; % O compensador deve contribuir com -phi para que o sistema seja compensado

%Parte para c�lculo do polo do compensador por avan�o de fase 

aux=1; % Variavel auxiliar para vetor de resultados
cont=1; % Variavel auxiliar para vetor do compensador
resultados=[]; % Declara��o do vetor de resultados para receber os valores relacionados ao atraso 
compensador=[]; % Declara��o do vetor de resultados para receber os valores relacionados ao compensador como um todo

for zav = -0.5: -0.1: -2 %Zero de avan�o escolhido na faixa de -0.5 a -2 com intervalos de 0.1
 
    angpz=angle(polosDes(1)-zav);% C�lculo do �ngulo entre o polo e zero
    beta_ang=angpz-Cont_angular; %beta_ang: �ngulo entre o polo de avan�o e polo desejado
    %Pela rela��o entre os �ngulos do tri�ngulo temos que 
    % Pavan�o= imag(polo desejado)/tan(angulo entre pAvan�o e polo
    % desejado)+abs(real(polo desejado))
    pav=-(wd/(abs(tan(beta_ang)))+sigma); % Como utilizei zav negativo, pav tamb�m � negativo
    
    % Pela condi��o de m�dulo temos que 
    % abs(Kc*Gavanco*G(s))=1
    % Determinamos ent�o o Kc
    Kc = 1/(abs((polosDes(1)-zav)/(polosDes(1)-pav))*abs(7/(polosDes(1)*(polosDes(1)+1.5)*(polosDes(1)+2.5))));
    Gav=Kc*(s-zav)/(s-pav); %Fun��o de transfer�ncia do compensador avan�o
    
    % Como s+1/T1=0 quando s=zav, isolando T1 temos:
    T1=1/(-zav); 
    % Como s+gama/T1=0 quando s=pav, isolando gama temos:
    gama=-pav*T1;
    
    %Projeto do compensador atraso
    %Kn=Kv/5
    %Kn=Kv*Kc/gama*beta
    beta=5*gama/Kc; %F�rmula do beta
    
    %Escolher T2 de modo que a condi��o de �ngulo e modulo sejam
    %satisfeitas
    
    for T2=1:1:10 %Varredura do T2
        anguloAtraso=-rad2deg(angle(((polosDes(1))+(1/T2))/((polosDes(1))+(1/(beta*T2))))); %C�lculo do �ngulo do compensador atraso
        moduloAtraso=abs((((polosDes(1))+(1/T2))/((polosDes(1))+(1/(beta*T2))))); %C�lculo do m�dulo do compensador atraso
     
        if anguloAtraso>0 && anguloAtraso<=5 && moduloAtraso>=0.96 %Condi��es de m�dulo e �ngulo
            Gatraso=(s+1/T2)/(s+1/(beta*T2)); %Fun��o de transfer�ncia do compensador atraso
            Gc=Gav*Gatraso; %Fun��o de transfer�ncia do compensador avan�o-atraso
            
            FTMAcompensado=FTMA*Gc; %Fun��o de transfer�ncia de malha aberta do sistema compensado
            FTMFcompensado=feedback(FTMAcompensado,1); %Fun��o de transfer�ncia de malha aberta do sistema compensado
            
            t=0:0.01:50;% Tempo da simula��o para entrada em degrau
            rdegrau=step(FTMFcompensado,t);%Com a FTMF compensada, aplica-se um degrau na entrada pela fun��o step
            
            m2=max(rdegrau); %Posi��o do m�ximo sobressinal da resposta ao degrau
            MP2=(m2-1)*100; %M�ximo sobressinal percentual
            
            if m2<1.36 %Condi��o do projeto para sobressinal de at� 36% acima
                
                %Busca de tr�s para frente para fazer menor esfor�o
                %computacional
                %Para isso, definimos o comprimento do vetor de tempo e
                %fazemos a busca de tr�s para frente
                comp=length(t); %fun��o length
                
                % Busca inversa pela faixa 5%
                while rdegrau(comp)>0.95 &&rdegrau(comp)<1.05
                    comp=comp-1;
                end
                %Atrav�s da posi��o de comp se calcula o tempo de
                %assentamento
                tassentamento=(comp-1)*0.01 ; %Posi��o no vetor multiplicado pelo tamanho dos intervalos
                if tassentamento <= 3.8 %Condi��o do tempo de estabiliza��o
                    resultados(aux,:)= [T2 m2 tassentamento]; %Armazenamento na matriz resultados dos valores de T2, m2 e tassentamento que satisfazem as restri��es
                    
                    aux=aux+1; %itera��o da vari�vel auxiliar
                end
            end
        end
    end
     if length(resultados)>0 % Verificar se houve algum resultado que satisfez as condi��es (importante porque ainda est� dentro de um loop)
        resultados_classificados = sortrows(resultados, 2); %Classificar os sobressinais(coluna 2) de modo crescente
        %Armazenamento na matriz compensador os valores relacionados ao
        %compensador que s�o op��es
        compensador(cont,:) = [angpz pav zav Kc T1 gama beta resultados_classificados(1,1) resultados_classificados(1,2) resultados_classificados(1,3)]; 
                      % Ordem [angpz pav zav Kc T1 gama beta T2 menorsobressinal tassentamento]
        cont=cont+ 1;
     end

end
%Agora definimos o melhor compensador avan�o e atraso depois de ter
%determinado o melhor atraso para cada avan�o 
%Como a matriz est� ordenada por menor sobressinal, foi escolhido a linha
% 1 da matriz como escolha de projeto por ter o melhor sobressinal

%Escolha dos par�metros do compensador
Kc_escolhido = compensador(1, 4);
T1_escolhido = compensador(1, 5);
T2_escolhido = compensador(1, 8);
gama_escolhido = compensador(1, 6);
beta_escolhido = compensador(1, 7);

Gav_escolhido=Kc_escolhido*((s+(1/T1_escolhido))/(s+(gama_escolhido/T1_escolhido))); %Fun��o de transfer�ncia do compensador avan�o escolhido
Gatraso_escolhido=((s+(1/T2_escolhido))/(s+(1/(beta_escolhido*T2_escolhido)))); %Fun��o de transfer�ncia do compensador atraso escolhido
Gc_escolhido=Gav_escolhido*Gatraso_escolhido;%Fun��o de transfer�ncia do compensador avan�o-atraso escolhido

FTMA_escolhido=FTMA*Gc_escolhido; %Fun��o de transfer�ncia de malha aberta do sistema compensado escolhido
FTMF_escolhido=feedback(FTMA_escolhido,1);%Fun��o de transfer�ncia de malha fechada do compensador escolhido

%Root Locus do sistema compensado com polos dominantes de ganho unit�rio
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
% Defini��o do t�tulo do gr�fico e do eixos
title('Grafico da resposta ao degrau dos sistema compensado');
xlabel('Tempo');
ylabel('Amplitude');

% Resposta a rampa unit�ria do sistema compensado com erro de regime
% permanente 5 vezes menor que o erro calculado e obtido no gr�fico do item
% D
figure(6)
hold on;
t=0:0.1:200;
c_rampa=lsim(FTMF_escolhido,t,t); % Aplica��o da rampa � FT de malha fechada do sistema n�o comoensado no intervalo de tempo t2
plot(t,c_rampa,'-');
plot(t,t);
grid;
% Defini��o do t�tulo do gr�fico e do eixos
title('Grafico da resposta a rampa unitaria dos sistema compensado');
xlabel('Tempo');
ylabel('Amplitude');
erro_rampaComp=(max(t)-max(c_rampa)); %Verificando que erro da rampa diminuiu 5 vezes
%Novo erro da entrada a rampa em regime permanente
%ess=0.107142857525417;