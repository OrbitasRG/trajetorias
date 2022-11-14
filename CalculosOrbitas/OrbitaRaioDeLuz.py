from http.client import HTTPResponse
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from scipy.integrate import quad
import math
from matplotlib.animation import FuncAnimation, writers
from matplotlib.patches import Circle

class OrbitaRaioDeLuz:
    # d = Escolha o valor do parâmetro de impacto $d$ (em km): 
    def calcular_raio_luz(self, v0):
        d = eval(v0)

        rs_sun = 3  # Hipótese sobre o raio de Schwarzschild do corpo central (em km). Ligeiramente maior que o do Sol.
        par_imp = 2.0 * d / rs_sun  # b
        k = 1 / (par_imp ** 2)

        rst = 50
        norbit = 10

        def w(u):
            w = u ** 2 - 2 * (u ** 3)
            return w

        umax = 1 / 3
        wmax = 1 / 27
        ust = 1 / rst

        coef = [-2, 1, 0, -k]
        roots = np.roots(coef)
        tp2 = roots[1]
        tp3 = roots[0]

        eps = 0.000000001

        if k < wmax and ust < umax:
            uint = ust
            uext = tp2 * (1 - eps)
            norbit = 1
        elif k > wmax:
            uint = ust
            uext = 0.5 * (1 - eps)
            norbit = 0.5
        elif k < wmax and ust > umax:
            uint = 0.5
            uext = tp3 * (1 + eps)
            norbit = 0.5
        else:
            print("Ha uma incoerencia entre os parametros fornecidos")

        v = sp.Symbol('v')

        def lambda_integrand(v):
            lambda_integrand = 1 / (v ** 2) * (k - w(v)) ** (-1 / 2)
            return lambda_integrand

        npoints = 300
        lambdatotal, errolambda = quad(lambda_integrand, uint, uext)  # Computes total affine parameter to go from uint to uext.
        dlambda = lambdatotal / npoints  # Sets the "time" step as 1/100 of the total "time".
        ud = [uext]
        for i in range(npoints + 20):
            ud.append(ud[i] - dlambda * ud[i] ** 2 * (k - w(ud[i])) ** (1 / 2.0))
            if ud[-1].imag != 0 or math.isnan(ud[-1]):
                ud = ud[:-2]
                break
        uc = ud[::-1]
        n = len(uc)


        def theta(v):
            theta = (k - w(v)) ** (-1 / 2)
            return theta


        delphi, erro = quad(theta, uint, uext)

        phi1 = []
        for i in range(len(uc)):
            a = quad(theta, uint, uc[i])
            phi1.append(abs(a[0]))

        phi2 = []
        for j in range(len(ud)):
            b = quad(theta, uext, ud[j])
            phi2.append(abs(b[0]))

        if norbit == 0.5:
            utotal = uc
        else:
            utotal = utotal = np.concatenate([uc, ud] * (norbit))

        accphi = [0] * (len(utotal))
        phi0 = np.arcsin(par_imp / rst)

        if norbit == 0.5:
            accphi = phi1
            x = [0] * (len(uc))
            y = [0] * (len(uc))
            for i in range(len(uc)):
                x[i] = (math.cos(phi0 + accphi[i])) / utotal[i] * (rs_sun / 2.0)
                y[i] = (math.sin(phi0 + accphi[i])) / utotal[i] * (rs_sun / 2.0)
        else:
            for i in range(norbit):
                for j in range(n):
                    accphi[j + (2 * i * n)] = 2 * i * delphi + phi1[j]
                    accphi[j + ((2 * i + 1) * n)] = ((2 * i) + 1) * delphi + phi2[j]
            x = [0] * (2 * norbit * n)
            y = [0] * (2 * norbit * n)
            for i in range(2 * norbit * n):
                x[i] = (math.cos(phi0 + accphi[i])) / utotal[i] * (rs_sun / 2.0)
                y[i] = (math.sin(phi0 + accphi[i])) / utotal[i] * (rs_sun / 2.0)

        fig = plt.figure()

        plt.xlabel("x (km)")
        plt.ylabel("y (km)")
        plt.gca().set_aspect('equal')
        ax = plt.gca()
        ax.spines['bottom'].set_color('white')
        ax.tick_params(axis='x', colors='white')
        ax.tick_params(axis='y', colors='white')
        ax.spines['top'].set_color('white')
        ax.spines['right'].set_color('white')
        ax.spines['left'].set_color('white')
        ax.xaxis.label.set_color('white')
        ax.yaxis.label.set_color('white')
        fig.patch.set_facecolor('#0E1117')
        ax.set_facecolor("black")
        circle = Circle((0, 0), rs_sun, color='dimgrey')
        plt.gca().add_patch(circle)
        plt.axis([- 0.8 * (rs_sun / 2.0) / uint, 0.8 * (rs_sun / 2.0) / uint, - 0.8 * (rs_sun / 2.0) / uint,
                0.8 * (rs_sun / 2.0) / uint])

        # Montagem do gif

        graph, = plt.plot([], [], 'k--', color="gold", markersize=3)

        def animate(i):
            graph.set_data(x[:i], y[:i])
            return graph,


        skipframes = int(len(x) / 200)
        if skipframes == 0:
            skipframes = 1

        ani2 = FuncAnimation(fig, animate, frames=range(0, len(x), skipframes), interval=10, blit=True,repeat=False)

        response = HTTPResponse(mimetype="image/png")
        plt.savefig(response, format="png")
        return response