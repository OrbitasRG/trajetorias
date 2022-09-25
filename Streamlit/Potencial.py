import streamlit as st
from streamlit_option_menu import option_menu
#import hydralit_components as hc
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from scipy.integrate import quad
import math
import time
from PIL import Image
from matplotlib.patches import Circle
import streamlit.components.v1 as components

st.set_page_config(
     page_title="Órbitas Relativisticas",
     page_icon=":comet:"
) 

hide_menu_style = """
        <style>
        #MainMenu {visibility: hidden;}
        </style>
        """
st.markdown(hide_menu_style, unsafe_allow_html=True)

selected = option_menu(menu_title=None, options=["Introdução","Corpos Massivos", "Raios de Luz"],icons = ["card-text","box","brightness-high"], default_index=0, orientation='horizontal',
)

#selected = option_menu(options=["Introdução", "Corpos Massivos", "Raios de Luz"], icon=["card-text","box","brightness-high"], default_index=0, orientation='horizontal',
#)


#st.set_page_config(page_title="Órbitas Relativísticas", page_icon=":star:")


#app = hy.HydraApp()


#@app.addapp(title='Introdução', icon="📜")
#def my_home():
if selected == "Introdução":
    st.title("Órbitas de corpos massivos")
    st.subheader('Caso newtoniano')
    st.write("De acordo com a lei da gravitação universal de Newton, o campo gravitacional externo a um corpo esférico de massa  $M$  (situado na origem do sistema de coordenadas) é")
    eq1 = r'''
    $$
    \vec{g} = - \frac{GM}{r^2} \hat{r} \qquad \qquad 
    $$
    '''
    st.write(eq1)
    st.write("onde $G = 6,67  •  10^{-11} m^3 kg^{-1} s^{-2}$ é a constante da gravitação universal e $r$ é a distância ao corpo central.")
    st.write("Uma partícula de massa $m$ sujeita a esse campo gravitacional descreve uma trajetória que é restrita a um plano e pode ser descrita pelas fuções $r(t)$ e $θ(t)$, que satisfazem")
    eq2 = r'''
    $$
    \frac{m}{2} \left(\frac{dr}{dt}\right)^2 + U_{efetiva} = E \qquad \qquad 
    $$
    '''
    st.write(eq2)
    st.write("e")
    eq3 = r'''
    $$
    \frac{d\theta}{dt} = \frac{L}{mr^2} \qquad \qquad 
    $$
    '''
    st.write(eq3)
    st.write("Aqui, $L$ é o módulo do momento angular, $E$ é a energia total e a energia potencial efetiva $U_{efetiva}$ é dada por")
    eq4 = r'''
    $$
    U_{efetiva} = -\frac{GMm}{r} + \frac{L^2}{2 m r^2} \qquad \qquad
    $$
    '''
    st.write(eq4)
    st.subheader("Caso Relativístico")
    st.write("A relatividade geral é a teoria atual que melhor descreve os fenômenos gravitacionais. Nessa teoria, a gravitação é descrita como resultado da curvatura do espaço-tempo. Uma das soluções mais simples, descrevendo a geometria do espaço-tempo externa a um corpo esférico, estático, de massa $M$ é a chamada solução de Schwarzschild. Partículas que orbitam o corpo central seguem trajetórias que também podem ser descritas pela equação (2), mas com a nova energia potencial efetiva:")
    eq5 = r'''
    $$
    U_{efetiva}^{(R)} = - \frac{GMm}{r} + \frac{L^2}{2 m r^2} - \frac{GML^2}{m c^2 r^3} \qquad \qquad 
    $$
    '''
    st.write(eq5) 

    st.title("Órbitas de raios de luz")
    st.write("Um raio de luz, no espaço-tempo de Schwarzschild, descreve uma trajetória que também pode ser escrita na forma de um problema unidimensional efetivo. Temos:")
    eq6 = r'''
    $$
    \frac{1}{\ell^2} \left(\frac{dr}{d\lambda}\right)^2 + V_{efetivo} = \frac{1}{d^2} \qquad \qquad
    $$
    '''
    st.write(eq6)
    st.write("com o potencial efetivo")
    eq7 = r'''
    $$
    V_{efetivo} = \frac{1}{r^2} \left( 1 - \frac{r_g}{r} \right) \qquad \qquad
    $$
    '''
    st.write(eq7)
    st.write("Aqui, $r_g = 2GM/c^2$ depende da massa do corpo central.")
    st.write("Nesse caso, o que cumpre o papel do parâmetro de 'energia' é o fator $1/d^2$, onde $d$ é o parâmetro de impacto do fóton, que pode ser obtido a partir da sua energia $e$ e momento angular $\ell$ como $d = c |\ell/e|$. Para entender o que d representa, considere um raio de luz que vem do infinito, se movendo paralelamente ao eixo $x$: o parâmetro de impacto $d$ é justamente a distância ao eixo-$x$, como mostra a figura abaixo.")
    image = Image.open('parametro.jpg')
    st.image(image)
    
    st.title("Próximos passos")
    st.write("Este simulador se destina ao estudo de órbitas relativísticas. Nele você pode simular o potencial efetivo e esboço de órbitas de corpos massivos e/ou raios de luz. Basta ir ao menu na parte superior desta página e escolher a partícula a ser analisada. Bom estudo!")
    

#@app.addapp(title='Corpos massivos', icon="🪨")
#def app2():
elif selected == "Corpos Massivos":
    st.title("Simulador para análise das órbitas relativísticas de corpos com massa")
    st.write("Funcionamento do programa: Seguindo os comandos abaixo, você deverá, primeiramente, inserir o valor do momento angular adimensional $Lc/(GMm)$. Em seguida, o programa exibirá o gráfico da energia potencial efetiva (adimensional): $U_{eff}^{(R)}/(mc^2)$. Você deverá então inserir o valor do parâmetro de energia da partícula teste, $E = E/(mc^2)$, que pode assumir qualquer valor maior que o mínimo de $U_{eff}^{(R)}$. O programa então retornará um gráfico correspondente à trajetória da partícula com esses parâmetros de energia e momento angular, para um corpo central ultracompacto com massa igual à do Sol.")
    st.subheader("Insira o valor do momento angular adimensional:")   
    momento = st.slider("Escolha entre 0 e 15",min_value=0.0, max_value=15.0, step = 0.1) 
    result1 = st.button("Gerar Potencial")

    if st.session_state.get('button') != True:
        st.session_state['button'] = result1

    if st.session_state['button'] == True:
        import numpy as np
        import matplotlib.pyplot as plt
        import sympy as sp
        from scipy.integrate import quad
        import math 
        from matplotlib.patches import Circle
     
        if momento==0:
               momento=0.000001
        else:
          momento=momento

        l = momento 

        def v(u, l):
            v = -u + (l ** 2) * (u ** 2) / 2 - (l ** 2) * (u ** 3)
            return v 

        if l > np.sqrt(12):
            coef = [0, - 3 * (l ** 2), l ** 2, -1]
            a = np.roots(coef)
            umax=sorted(a)[1].real
            umin=sorted(a)[0].real
            vmin = v(umin, l)
            vmax = v(umax, l)
            vlim = vmax
            rmax = 2/umin
            st.write("O mínimo e o máximo da energia potencial efetiva são",vmin,"e",vmax," (ver pontos no gráfico)")
            st.write("")
            r = np.arange(2, rmax, rmax/30000)
            u = 1 / r
            fig1 = plt.figure()
            plt.subplot(1, 1, 1)
            plt.axhline(0, linewidth=0.3, color='white')
            plt.plot(1.477*r, v(u, l), color="white")
            plt.plot([1.477/umin,1.477/umax],[vmin,vmax],'bo', color="gold")
            plt.xlabel("r [km]")
            plt.axis([-1, 1.2*rmax, -0.5, vlim + 0.1])
            #plt.axis([1.477/umax*0.05, 1.477/umin*1.2, -0.5, vlim*1.2])
            ax = plt.gca()
            ax.spines['bottom'].set_color('white')
            ax.tick_params(axis='x', colors='white')
            ax.tick_params(axis='y', colors='white')
            ax.spines['top'].set_color('white') 
            ax.spines['right'].set_color('white')
            ax.spines['left'].set_color('white')
            ax.xaxis.label.set_color('white')
            ax.yaxis.label.set_color('white')
            fig1.patch.set_facecolor('#0E1117')
            ax.set_facecolor("black")
            plt.show()
            st.pyplot(fig1)
            
            st.subheader("Insira o valor do parâmetro de energia:") 
            E = st.number_input('Escolha um valor maior que o mínimo da energia potencial efetiva; calculada no passo anteior', step=1e-5, format="%.5f")
            result111 = st.button("Continuar")
            
            if st.session_state.get('button11') != True:
                       st.session_state['button11'] = result111
                       
            if st.session_state['button11'] == True:           
               if E<vmin:
                    st.write("ATENÇÃO: Você escolheu um valor de parâmetro de energia menor que o valor mínimo permitido. Escolha novamente um valor dentro do limite especificado antes de continuar.")   
               else:
                   st.subheader("Insira o número de voltas completas na órbita:")   
                   st.write("Para uma órbita ligada ($U_{eff,min} ≤ E < 0$), escolha também o número de órbitas que deseja traçar.")
                   norbit = st.slider("Escolha entre 1 e 20",min_value=1, max_value=20, step = 1)
                   result11 = st.button("Gerar Órbita")

                   if st.session_state.get('button1') != True:
                            st.session_state['button1'] = result11

                   if st.session_state['button1'] == True:
                           import numpy as np
                           import matplotlib.pyplot as plt
                           import sympy as sp
                           from scipy.integrate import quad
                           import math 
                           from matplotlib.patches import Circle
                           if E==0:
                               E=E+1e-10    
                           coef = [- (l ** 2), (l ** 2) / 2, -1, -E]
                           roots = np.roots(coef)
                           tp1 = roots[2]
                           tp2 = roots[1]
                           tp3 = roots[0]

                           eps = 0.00000001
                           rst = 10*l/(E+0.5)
                           ust = 1 / rst
                           correction = 1.477

                           if l > math.sqrt(12):
                               if E < 0 and ust < tp2.real:
                                   u1 = tp1.real * (1 + eps)
                                   u2 = tp2.real * (1 - eps)
                               elif 0 < E < vmax and ust < tp2.real:
                                   u1 = ust
                                   u2 = tp2.real * (1 - eps)
                                   norbit = 1
                               elif E < vmax and ust > tp3.real:
                                   u1 = 0.5
                                   u2 = tp3.real * (1 + eps)
                                   norbit = 0.5
                               elif E > vmax:
                                   u1 = ust
                                   u2 = 0.5
                                   norbit = 0.5
                           else:
                               if E >= 0:
                                   u1 = ust
                                   u2 = 0.5
                                   norbit = 0.5
                               else:
                                   u1 = tp1.real * (1 + eps)
                                   u2 = 0.5
                                   norbit = 0.5

                           def theta(w):
                               theta = (l/(2**(1/2)))*((E-v(w,l))**(-1/2))
                               return theta

                           delphi, erro = quad(theta, u1, u2)

                           n = 1000
                           #ui=np.arange(u1,u2/2,(u2/2 - u1)/n)
                           #uf=np.arange(u2/2, u2, (u2 - u2/2)/(10*n))
                           #uc=ui.extend(uf)
                           #ud=reverse(uc)
                           uc = np.arange(u1, u2, (u2 - u1)/n)
                           ud = np.arange(u2, u1, (u1 - u2)/n)

                           phi1 = []
                           for i in range(len(uc)):
                               a = quad(theta, u1, uc[i])
                               phi1.append(abs(a[0]))

                           phi2 = []
                           for j in range(len(ud)):
                               b = quad(theta, u2, ud[j])
                               phi2.append(abs(b[0]))

                           if norbit == 0.5:
                               utotal = uc
                           else:
                               utotal = np.concatenate([uc, ud]*(norbit))

                           accphi = [0]*(len(utotal))

                           if norbit == 0.5:
                               accphi = phi1
                               x = [0] * (len(uc))
                               y = [0] * (len(uc))
                               for i in range(len(uc)):
                                   x[i] = (math.cos(accphi[i])) / utotal[i] *correction
                                   y[i] = (math.sin(accphi[i])) / utotal[i] *correction
                           else:
                               for i in range (norbit):
                                   for j in range (n):
                                       accphi[j+(2*i*n)] = 2 * i * delphi + phi1[j]
                                       accphi[j+((2*i+1)*n)] = ((2*i)+1)*delphi + phi2[j]
                               x = [0] * (2 * norbit * n)
                               y = [0] * (2 * norbit * n)
                               for i in range(2 * norbit * n):
                                   x[i] = ((math.cos(accphi[i])) / utotal[i])*correction
                                   y[i] = ((math.sin(accphi[i])) / utotal[i])*correction

                           fig2 = plt.figure() 
                           plt.plot(x, y, color="gold")
                           plt.xlabel("x [km]")
                           plt.ylabel("y [km]")
                           circle = Circle((0,0), 2*1.477, color = 'dimgrey')
                           plt.gca().add_patch(circle)
                           plt.gca().set_aspect('equal')
                           plt.axis([(-1 / u1-1)*1.477 , (1 / u1+1)*1.477, (-1 / u1-1)*1.477 , (1/u1 +1)*1.477])
                           ax = plt.gca()
                           ax.spines['bottom'].set_color('white')
                           ax.tick_params(axis='x', colors='white')
                           ax.tick_params(axis='y', colors='white')
                           ax.spines['top'].set_color('white') 
                           ax.spines['right'].set_color('white')
                           ax.spines['left'].set_color('white')
                           ax.xaxis.label.set_color('white')
                           ax.yaxis.label.set_color('white')
                           fig2.patch.set_facecolor('#0E1117')
                           ax.set_facecolor("black")
                           plt.show()
                           st.pyplot(fig2)

                           st.session_state['button1'] = False                                                                             
            

        else:
                 vlim = 0
                 rmax = 80
                 st.write("Para esse valor do momento angular, a energia potencial efetiva não possui um mínimo ou máximo local.")
                 st.write("")
                 r = np.arange(2, rmax, rmax/300)
                 u = 1 / r
                 fig1 = plt.figure()
                 plt.subplot(1, 1, 1)
                 plt.plot(1.477*r, v(u, l), color="white")
                 plt.axhline(0, linewidth=0.3, color='white')
                 plt.xlabel("r [km]")
                 plt.axis([0, 1.477*rmax, -0.5, vlim + 0.1])
                 ax = plt.gca()
                 ax.spines['bottom'].set_color('white')
                 ax.tick_params(axis='x', colors='white')
                 ax.tick_params(axis='y', colors='white')
                 ax.spines['top'].set_color('white') 
                 ax.spines['right'].set_color('white')
                 ax.spines['left'].set_color('white')
                 ax.xaxis.label.set_color('white')
                 ax.yaxis.label.set_color('white')
                 fig1.patch.set_facecolor('#0E1117')
                 ax.set_facecolor("black")
                 plt.show()
                 st.pyplot(fig1)

                 st.subheader("Insira o valor do parâmetro de energia:") 
                 E = st.number_input('Escolha um valor maior que o mínimo da energia potencial efetiva; calculada no passo anteior', step=1e-5, format="%.5f")
                 result111 = st.button("Continuar")

                 if st.session_state.get('button11') != True:
                            st.session_state['button11'] = result111

                 if st.session_state['button11'] == True:
                             st.subheader("Insira o número de voltas completas na órbita:")   
                             st.write("Para uma órbita ligada ($U_{eff,min} ≤ E < 0$), escolha também o número de órbitas que deseja traçar.")
                             norbit = st.slider("Escolha entre 1 e 20",min_value=1, max_value=20, step = 1)
                             result11 = st.button("Gerar Órbita")

                             if st.session_state.get('button1') != True:
                                      st.session_state['button1'] = result11

                             if st.session_state['button1'] == True:
                                     import numpy as np
                                     import matplotlib.pyplot as plt
                                     import sympy as sp
                                     from scipy.integrate import quad
                                     import math 
                                     from matplotlib.patches import Circle
                                     if E==0:
                                         E=E+1e-10    
                                     coef = [- (l ** 2), (l ** 2) / 2, -1, -E]
                                     roots = np.roots(coef)
                                     tp1 = roots[2]
                                     tp2 = roots[1]
                                     tp3 = roots[0]

                                     eps = 0.00000001
                                     rst = 10*l/(E+0.5)
                                     ust = 1 / rst
                                     correction = 1.477

                                     if l > math.sqrt(12):
                                         if E < 0 and ust < tp2.real:
                                             u1 = tp1.real * (1 + eps)
                                             u2 = tp2.real * (1 - eps)
                                         elif 0 < E < vmax and ust < tp2.real:
                                             u1 = ust
                                             u2 = tp2.real * (1 - eps)
                                             norbit = 1
                                         elif E < vmax and ust > tp3.real:
                                             u1 = 0.5
                                             u2 = tp3.real * (1 + eps)
                                             norbit = 0.5
                                         elif E > vmax:
                                             u1 = ust
                                             u2 = 0.5
                                             norbit = 0.5
                                     else:
                                         if E >= 0:
                                             u1 = ust
                                             u2 = 0.5
                                             norbit = 0.5
                                         else:
                                             u1 = tp1.real * (1 + eps)
                                             u2 = 0.5
                                             norbit = 0.5

                                     def theta(w):
                                         theta = (l/(2**(1/2)))*((E-v(w,l))**(-1/2))
                                         return theta

                                     delphi, erro = quad(theta, u1, u2)

                                     n = 1000
                                     #ui=np.arange(u1,u2/2,(u2/2 - u1)/n)
                                     #uf=np.arange(u2/2, u2, (u2 - u2/2)/(10*n))
                                     #uc=ui.extend(uf)
                                     #ud=reverse(uc)
                                     uc = np.arange(u1, u2, (u2 - u1)/n)
                                     ud = np.arange(u2, u1, (u1 - u2)/n)

                                     phi1 = []
                                     for i in range(len(uc)):
                                         a = quad(theta, u1, uc[i])
                                         phi1.append(abs(a[0]))

                                     phi2 = []
                                     for j in range(len(ud)):
                                         b = quad(theta, u2, ud[j])
                                         phi2.append(abs(b[0]))

                                     if norbit == 0.5:
                                         utotal = uc
                                     else:
                                         utotal = np.concatenate([uc, ud]*(norbit))

                                     accphi = [0]*(len(utotal))

                                     if norbit == 0.5:
                                         accphi = phi1
                                         x = [0] * (len(uc))
                                         y = [0] * (len(uc))
                                         for i in range(len(uc)):
                                             x[i] = (math.cos(accphi[i])) / utotal[i] *correction
                                             y[i] = (math.sin(accphi[i])) / utotal[i] *correction
                                     else:
                                         for i in range (norbit):
                                             for j in range (n):
                                                 accphi[j+(2*i*n)] = 2 * i * delphi + phi1[j]
                                                 accphi[j+((2*i+1)*n)] = ((2*i)+1)*delphi + phi2[j]
                                         x = [0] * (2 * norbit * n)
                                         y = [0] * (2 * norbit * n)
                                         for i in range(2 * norbit * n):
                                             x[i] = ((math.cos(accphi[i])) / utotal[i])*correction
                                             y[i] = ((math.sin(accphi[i])) / utotal[i])*correction

                                     fig2 = plt.figure() 
                                     plt.plot(x, y, color="gold")
                                     plt.xlabel("x [km]")
                                     plt.ylabel("y [km]")
                                     circle = Circle((0,0), 2*1.477, color = 'dimgrey')
                                     plt.gca().add_patch(circle)
                                     plt.gca().set_aspect('equal')
                                     plt.axis([(-1 / u1-1)*1.477 , (1 / u1+1)*1.477, (-1 / u1-1)*1.477 , (1/u1 +1)*1.477])
                                     ax = plt.gca()
                                     ax.spines['bottom'].set_color('white')
                                     ax.tick_params(axis='x', colors='white')
                                     ax.tick_params(axis='y', colors='white')
                                     ax.spines['top'].set_color('white') 
                                     ax.spines['right'].set_color('white')
                                     ax.spines['left'].set_color('white')
                                     ax.xaxis.label.set_color('white')
                                     ax.yaxis.label.set_color('white')
                                     fig2.patch.set_facecolor('#0E1117')
                                     ax.set_facecolor("black")
                                     plt.show()
                                     st.pyplot(fig2)

                                     st.session_state['button1'] = False

#@app.addapp(title='Raios de luz', icon="💡")
#def app2():
elif selected == "Raios de Luz":
    st.title("Simulador para análise das órbitas relativísticas de raios de luz")
    st.write("Desta vez, a forma do potencial é fixa e está ilustrada abaixo. O máximo do potencial efetivo ocorre para $r = 1.5 r_g$.")
    imagepot = Image.open('potencialluz.png')
    st.image(imagepot)
    st.write("Funcionamento do programa: Seguindo os comandos abaixo, você deverá inserir o valor do parâmetro de impacto (em unidades do raio gravitacional do corpo central). Você também pode escolher $d≈2.5981$, o valor do parâmetro de impacto crítico ($d=(V_{max})^{-1/2}$). O programa então retornará um gráfico correspondente à trajetória da partícula com esse parâmetro de impacto. Em (1) você pode explorar uma faixa definida de valores de  $d$ ; em (2), você pode inserir o valor desejado diretamente.")
    st.subheader("(1) Escolha o valor do parâmetro de impacto $d$, em unidades de $r_g$:")   
    d = st.slider("Escolha entre 0 e 20",min_value=0.0, max_value=20.0, step = 0.1) 
    result2 = st.button("(1) Gerar Órbita")

    if st.session_state.get('button') != True:
        st.session_state['button'] = result2

    if st.session_state['button'] == True:
        import numpy as np
        import matplotlib.pyplot as plt
        import sympy as sp
        from scipy.integrate import quad
        import math
        from matplotlib.patches import Circle
     
        if d==0:
          d=0.00001
        else:
          d=d

        b = 2*d

        rst = 50
        norbit = 10

        r = np.arange(0.1, 80, 0.1)
        u = 1 / r
        k = 1/(b**2)

        def w(u):
            w = u**2 - 2*(u**3)
            return w

        umax = 1/3
        wmax = 1/27
        ust = 1/rst

        coef = [-2, 1 , 0, -k]
        roots = np.roots(coef)
        tp2 = roots[1]
        tp3 = roots[0]

        eps = 0.000000001

        if k == wmax:
            phi = np.arange(0, 2 * math.pi, math.pi / 500)
            rmax = [1 / umax] * len(phi)
            xymax = [[0] * (len(phi)), [0] * (len(phi))]
            xymin = [[0] * (len(phi)), [0] * (len(phi))]
            for i in range(len(phi)):
                xymax[0][i] = rmax[i] * math.cos(phi[i]) /2 #Divisão por 2 para colocar em unidades de rg, não de M.
                xymax[1][i] = rmax[i] * math.sin(phi[i]) /2
            plt.figure()
            plt.subplot(1, 1, 1)
            plt.xlabel("x / rg")
            plt.ylabel("y / rg")
            circle = Circle((0,0), 2*1.477, color = 'black')
            plt.gca().add_patch(circle)
            plt.gca().set_aspect('equal')
            plt.plot(xymax[0], xymax[1], 'k--', color='black')
            plt.axis([1 / umax + 1, -1 / umax - 1, 1 / umax + 1, -1 / umax - 1])
        else:
            if k < wmax and ust < umax:
                uint = ust
                uext = tp2*(1 - eps)
                norbit = 1
            elif k > wmax:
                uint = ust
                uext = 0.5 * (1 - eps)
                norbit=0.5
            elif k < wmax and ust > umax:
                uint = 0.5
                uext = tp3 * (1 + eps)
                norbit =  0.5
            else:
                print("Ha uma incoerencia entre os parametros fornecidos")

            v = sp.Symbol('v')

            def theta(v):
                theta = (k-w(v))**(-1/2)
                return theta

            delphi, erro = quad(theta, uint, uext)

            n = 500
            uc = np.arange(uint, uext, (uext - uint) / n)
            ud = np.arange(uext, uint, (uint - uext) / n)

            phi1 = []
            for i in range(len(uc)):
                a = quad(theta, uint, uc[i])
                phi1.append(abs(a[0]))

            phi2 = []
            for j in range(len(ud)):
                b = quad(theta, uext, ud[j])
                phi2.append(abs(b[0]))

            if norbit==0.5:
                utotal = uc
            else:
                utotal = utotal = np.concatenate([uc, ud]*(norbit))

            accphi = [0]*(len(utotal))

            if norbit == 0.5:
                accphi = phi1
                x = [0] * (len(uc))
                y = [0] * (len(uc))
                for i in range(len(uc)):
                    x[i] = (math.cos(accphi[i])) / utotal[i] /2 #Divisão por 2 para colocar em unidades de rg, não de M.
                    y[i] = (math.sin(accphi[i])) / utotal[i] /2
            else:
                for i in range(norbit):
                    for j in range(n):
                        accphi[j + (2 * i * n)] = 2 * i * delphi + phi1[j]
                        accphi[j + ((2 * i + 1) * n)] = ((2 * i) + 1) * delphi + phi2[j]
                x = [0] * (2 * norbit * n)
                y = [0] * (2 * norbit * n)
                for i in range(2 * norbit * n):
                    x[i] = (math.cos(accphi[i])) / utotal[i] /2 
                    y[i] = (math.sin(accphi[i])) / utotal[i] /2 #Divisão por 2 para colocar em unidades de rg, não de M.

            fig2 = plt.figure()
            plt.subplot(1, 1, 1)
            plt.plot(x, y, 'k--', color='gold')
            plt.xlabel("x [rg]")
            plt.ylabel("y [rg]")
            circle = Circle((0,0), 1, color = 'dimgrey')
            plt.gca().add_patch(circle)
            plt.gca().set_aspect('equal')
            plt.axis([(-1 / uint + 10)/2 , (1 / uint -10)/2 , (-1 / uint +10)/2 , (1 / uint -10)/2 ])
            ax = plt.gca()
            ax.spines['bottom'].set_color('white')
            ax.tick_params(axis='x', colors='white')
            ax.tick_params(axis='y', colors='white')
            ax.spines['top'].set_color('white') 
            ax.spines['right'].set_color('white')
            ax.spines['left'].set_color('white')
            ax.xaxis.label.set_color('white')
            ax.yaxis.label.set_color('white')
            fig2.patch.set_facecolor('#0E1117')
            ax.set_facecolor("black")
            plt.show()
            st.pyplot(fig2)
            
            st.session_state['button'] = False  
    
    st.subheader("(2) Escolha o valor do parâmetro de impacto $d$, em unidades de $r_g$:")   
    d2 = st.number_input('Insira um valor de parâmetro de impacto',step=1e-4, format="%.4f")
    result3 = st.button("(2) Gerar Órbita")

    if st.session_state.get('button2') != True:
        st.session_state['button2'] = result3

    if st.session_state['button2'] == True:
             import numpy as np
             import matplotlib.pyplot as plt
             import sympy as sp
             from scipy.integrate import quad
             import math
             from matplotlib.patches import Circle

             if d2==0:
               d2=0.00001
             else:
               d2=d2

             b = 2*d2

             rst = 50
             norbit = 10

             r = np.arange(0.1, 80, 0.1)
             u = 1 / r
             k = 1/(b**2)

             def w(u):
                 w = u**2 - 2*(u**3)
                 return w

             umax = 1/3
             wmax = 1/27
             ust = 1/rst

             coef = [-2, 1 , 0, -k]
             roots = np.roots(coef)
             tp2 = roots[1]
             tp3 = roots[0]

             eps = 0.000000001

             if k == wmax:
                 phi = np.arange(0, 2 * math.pi, math.pi / 500)
                 rmax = [1 / umax] * len(phi)
                 xymax = [[0] * (len(phi)), [0] * (len(phi))]
                 xymin = [[0] * (len(phi)), [0] * (len(phi))]
                 for i in range(len(phi)):
                     xymax[0][i] = rmax[i] * math.cos(phi[i]) /2 #Divisão por 2 para colocar em unidades de rg, não de M.
                     xymax[1][i] = rmax[i] * math.sin(phi[i]) /2
                 plt.figure()
                 plt.subplot(1, 1, 1)
                 plt.xlabel("x / rg")
                 plt.ylabel("y / rg")
                 circle = Circle((0,0), 2*1.477, color = 'black')
                 plt.gca().add_patch(circle)
                 plt.gca().set_aspect('equal')
                 plt.plot(xymax[0], xymax[1], 'k--', color='black')
                 plt.axis([1 / umax + 1, -1 / umax - 1, 1 / umax + 1, -1 / umax - 1])
             else:
                 if k < wmax and ust < umax:
                     uint = ust
                     uext = tp2*(1 - eps)
                     norbit = 1
                 elif k > wmax:
                     uint = ust
                     uext = 0.5 * (1 - eps)
                     norbit=0.5
                 elif k < wmax and ust > umax:
                     uint = 0.5
                     uext = tp3 * (1 + eps)
                     norbit =  0.5
                 else:
                     print("Ha uma incoerencia entre os parametros fornecidos")

                 v = sp.Symbol('v')

                 def theta(v):
                     theta = (k-w(v))**(-1/2)
                     return theta

                 delphi, erro = quad(theta, uint, uext)

                 n = 500
                 uc = np.arange(uint, uext, (uext - uint) / n)
                 ud = np.arange(uext, uint, (uint - uext) / n)

                 phi1 = []
                 for i in range(len(uc)):
                     a = quad(theta, uint, uc[i])
                     phi1.append(abs(a[0]))

                 phi2 = []
                 for j in range(len(ud)):
                     b = quad(theta, uext, ud[j])
                     phi2.append(abs(b[0]))

                 if norbit==0.5:
                     utotal = uc
                 else:
                     utotal = utotal = np.concatenate([uc, ud]*(norbit))

                 accphi = [0]*(len(utotal))

                 if norbit == 0.5:
                     accphi = phi1
                     x = [0] * (len(uc))
                     y = [0] * (len(uc))
                     for i in range(len(uc)):
                         x[i] = (math.cos(accphi[i])) / utotal[i] /2 #Divisão por 2 para colocar em unidades de rg, não de M.
                         y[i] = (math.sin(accphi[i])) / utotal[i] /2
                 else:
                     for i in range(norbit):
                         for j in range(n):
                             accphi[j + (2 * i * n)] = 2 * i * delphi + phi1[j]
                             accphi[j + ((2 * i + 1) * n)] = ((2 * i) + 1) * delphi + phi2[j]
                     x = [0] * (2 * norbit * n)
                     y = [0] * (2 * norbit * n)
                     for i in range(2 * norbit * n):
                         x[i] = (math.cos(accphi[i])) / utotal[i] /2 
                         y[i] = (math.sin(accphi[i])) / utotal[i] /2 #Divisão por 2 para colocar em unidades de rg, não de M.

                 fig3 = plt.figure()
                 plt.subplot(1, 1, 1)
                 plt.plot(x, y, 'k--', color='gold')
                 plt.xlabel("x [rg]")
                 plt.ylabel("y [rg]")
                 circle = Circle((0,0), 1, color = 'dimgrey')
                 plt.gca().add_patch(circle)
                 plt.gca().set_aspect('equal')
                 plt.axis([(-1 / uint + 10)/2 , (1 / uint -10)/2 , (-1 / uint +10)/2 , (1 / uint -10)/2 ])
                 ax = plt.gca()
                 ax.spines['bottom'].set_color('white')
                 ax.tick_params(axis='x', colors='white')
                 ax.tick_params(axis='y', colors='white')
                 ax.spines['top'].set_color('white') 
                 ax.spines['right'].set_color('white')
                 ax.spines['left'].set_color('white')
                 ax.xaxis.label.set_color('white')
                 ax.yaxis.label.set_color('white')
                 fig3.patch.set_facecolor('#0E1117')
                 ax.set_facecolor("black")
                 plt.show()
                 st.pyplot(fig3)

                 st.session_state['button2'] = False             

#Run the whole lot, we get navbar, state management and app isolation, all with this tiny amount of work.
#app.run()
