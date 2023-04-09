# Trabalho Computacional 2022.3 CEL038

Considere o sistema de controle abaixo.Para este sistema, faça com MatLab ou Python um código para:

$G(s)=\frac{7}{s(s+1.5)(s+2.5)}$

Item a) Construir o gráfico do Lugar das Raízes (LR). Destacar no gráfico do LR os polos de malha fechada de ganho
unitário e incluir a reta de amortecimento $\zeta$ e círculo de $\omega_n$ destes polos dominantes. Adicionar também ao LR
a reta de amortecimento 0,4 e o círculo de $\omega_n$ de 2 rad/s.

Item b) Aplicar um degrau unitário na referência r(t) e gerar o gráfico da resposta no tempo para a saída c(t), com 50s
de simulação. Utilize um vetor de tempo t, com intervalos de 0,01s , até 50s.

Item c) A partir da saída c(t) para o degrau unitário, obter, através do código, o valor do sobressinal percentual e da
frequência $\it{f}$ , com o valor $\omega_d$ correspondente. Compare o valor de ωd obtido com a parte imaginária dos polos
dominantes do item (a), usando os comentários para falar da relação entre ambos.

Item d) Aplicar uma rampa na referência r(t) e gerar o gráfico da resposta no tempo para a saída c(t), com 200s de
simulação. Calcular com o código o valor da constante de velocidade K$_v$ para G(s) e o erro de regime
permanente para a constante K$_v$ calculada. Compare o erro de regime permanente calculado com o erro
apresentado pelo gráfico da saída e comente a comparação.

Item e)Projetar, através do código, um compensador Gc(s) com $\gamma \neq \beta(\gamma > 1 \space e \space \beta > 1)$, em série com G(s), para
aumentar o amortecimento $\zeta$ dos polos dominantes para 0,4 e para ajustar a frequência natural não amortecida
$\omega_n$ para 2 rad/s. Este compensador também deve reduzir em 5 vezes o erro de regime permanente para entrada
em rampa. Fazer uma busca para a sintonia do compensador Gc(s), de forma que, para a entrada ao degrau
unitário, o sobressinal fique em até 36% e o tempo de assentamento (faixa de 5%) seja $\leq$ 3,8 s. Após projetar e
sintonizar Gc(s), traçar o LR do sistema compensado, destacando os polos de malha fechada dominantes que
atendem $\zeta$ e $\omega_n$ especificados. Traçar o gráfico da resposta ao degrau unitário em r(t) do sistema compensado
até 50s, demonstrando o sobressinal e o tempo de assentamento especificados. Traçar o gráfico da resposta a
rampa até 200s, demonstrando que o erro de velocidade reduziu 5 vezes após a compensação.

