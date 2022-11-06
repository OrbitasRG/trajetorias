import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from scipy.integrate import quad
import math
from matplotlib.patches import Circle

def w(u):
    return u ** 2 - 2 * (u ** 3)

class OrbitaPotencialRaioDeLuz:

    def __init__(self, impacto):
        d = eval(impacto)
        if d == 0:
            d = 0.00001

        self.b = 2 * d
        self.k = 1 / (self.b ** 2)
        self.rst = 50
        self.norbit = 10
        self.r = np.arange(0.1, 80, 0.1)
        self.u = 1 / self.r
        self.umax = 1 / 3
        self.wmax = 1 / 27
        self.ust = 1 / self.rst
        self.coef = [-2, 1, 0, -self.k]
        self.roots = np.roots(self.coef)
        self.tp2 = self.roots[1]
        self.tp3 = self.roots[0]
        self.eps = 0.000000001

    def calcularOrbita(self):
        if self.k == self.wmax:
            phi = np.arange(0, 2 * math.pi, math.pi / 500)
            rmax = [1 / self.umax] * len(phi)
            xymax = [[0] * (len(phi)), [0] * (len(phi))]
            xymin = [[0] * (len(phi)), [0] * (len(phi))]
            for i in range(len(phi)):
                xymax[0][i] = rmax[i] * math.cos(phi[i]) / 2  # Divisão por 2 para colocar em unidades de rg, não de M.
                xymax[1][i] = rmax[i] * math.sin(phi[i]) / 2

            #PLOT ORBITA 1

            plt.figure()
            plt.subplot(1, 1, 1)
            plt.xlabel("x / rg")
            plt.ylabel("y / rg")
            circle = Circle((0, 0), 2 * 1.477, color='black')
            plt.gca().add_patch(circle)
            plt.gca().set_aspect('equal')
            plt.plot(xymax[0], xymax[1], 'k--', color='black')
            plt.axis([1 / self.umax + 1, -1 / self.umax - 1, 1 / self.umax + 1, -1 / self.umax - 1])

        else:
            if self.k < self.wmax and self.ust < self.umax:
                uint = self.ust
                uext = self.tp2 * (1 - self.eps)
                norbit = 1
            elif self.k > self.wmax:
                uint = self.ust
                uext = 0.5 * (1 - self.eps)
                norbit = 0.5
            elif self.k < self.wmax and self.ust > self.umax:
                uint = 0.5
                uext = self.tp3 * (1 + self.eps)
                norbit = 0.5
            else:
                print("Ha uma incoerencia entre os parametros fornecidos")

            v = sp.Symbol('v')


            def theta(v):
                theta = (self.k - w(v)) ** (-1 / 2)
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

            if norbit == 0.5:
                utotal = uc
            else:
                utotal = utotal = np.concatenate([uc, ud] * (norbit))

            accphi = [0] * (len(utotal))

            if norbit == 0.5:
                accphi = phi1
                x = [0] * (len(uc))
                y = [0] * (len(uc))
                for i in range(len(uc)):
                    x[i] = (math.cos(accphi[i])) / utotal[
                        i] / 2  # Divisão por 2 para colocar em unidades de rg, não de M.
                    y[i] = (math.sin(accphi[i])) / utotal[i] / 2
            else:
                for i in range(norbit):
                    for j in range(n):
                        accphi[j + (2 * i * n)] = 2 * i * delphi + phi1[j]
                        accphi[j + ((2 * i + 1) * n)] = ((2 * i) + 1) * delphi + phi2[j]
                x = [0] * (2 * norbit * n)
                y = [0] * (2 * norbit * n)
                for i in range(2 * norbit * n):
                    x[i] = (math.cos(accphi[i])) / utotal[i] / 2
                    y[i] = (math.sin(accphi[i])) / utotal[
                        i] / 2  # Divisão por 2 para colocar em unidades de rg, não de M.

            #PLOT ORBITA 2

            fig2 = plt.figure()
            plt.subplot(1, 1, 1)
            plt.plot(x, y, 'k--', color='gold')
            plt.xlabel("x [rg]")
            plt.ylabel("y [rg]")
            circle = Circle((0, 0), 1, color='dimgrey')
            plt.gca().add_patch(circle)
            plt.gca().set_aspect('equal')
            plt.axis([(-1 / uint + 10) / 2, (1 / uint - 10) / 2, (-1 / uint + 10) / 2, (1 / uint - 10) / 2])
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
        


