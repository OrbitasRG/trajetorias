import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from scipy.integrate import quad
import math
from matplotlib.patches import Circle

def v(u, l):
    v = -u + (l ** 2) * (u ** 2) / 2 - (l ** 2) * (u ** 3)
    return v

class OrbitaCorposMassivos:
    umax = 0
    umin = 0
    vmin = 0
    vmax = 0
    vlim = 0
    rmax = 0
    a = 0
    r = 0
    u = 0
    l = 0.0
    root = []
    tp1 = 0
    tp2 = 0
    tp3 = 0

    def __init__(self, valor):
        momento = eval(valor)

        if momento == 0:
            momento = 0.000001
        else:
            momento = momento
        
        self.l = momento

        self.coef = [0, - 3 * (self.l ** 2), self.l ** 2, -1]
        self.a = np.roots(self.coef)
        self.umax = sorted(self.a)[1].real
        self.umin = sorted(self.a)[0].real
        self.vmin = v(self.umin, self.l)
        self.vmax = v(self.umax, self.l)
        self.vlim = self.vmax
        self.rmax = 2 / self.umin
        self.r = np.arange(2, self.rmax, self.rmax / 30000)
        self.u = 1 / self.r
   
    def calcularCorposMassivos(self):
        return self.l > np.sqrt(12)

    def calcularPlotPotencialVisaoUm(self):
        #PLOT DO POTENCIAL

        fig1 = plt.figure()
        plt.subplot(1, 1, 1)
        plt.axhline(0, linewidth=0.3, color='white')
        plt.plot(1.477 * self.r, v(self.u, self.l), color="white")
        plt.plot([1.477 / self.umin, 1.477 / self.umax], [self.vmin, self.vmax], 'bo', color="gold")
        plt.xlabel("r [km]")
        plt.axis([-1, 1.2 * self.rmax, -0.5, self.vlim + 0.1])
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
        return plt

    def calcularPlotPotencialVisaoDois(self, energia, numeroOrbitas):
        E = eval(energia)
        norbit = eval(numeroOrbitas)

        if norbit < 1 or norbit > 20:
            print('Escala deve estar entre 1 e 20')
            return

        if E < self.vmin:
            return
       
        if E == 0:
            E = E + 1e-10
        
        self.coef = [- (self.l ** 2), (self.l ** 2) / 2, -1, -E]
        self.roots = np.roots(self.coef)
        self.tp1 = self.roots[2]
        self.tp2 = self.roots[1]
        self.tp3 = self.roots[0]

        eps = 0.00000001
        rst = 10 * self.l / (E + 0.5)
        ust = 1 / rst
        correction = 1.477

        if self.l > math.sqrt(12):
            if E < 0 and ust < self.tp2.real:
                u1 = self.tp1.real * (1 + eps)
                u2 = self.tp2.real * (1 - eps)
            elif 0 < E < self.vmax and ust < self.tp2.real:
                u1 = ust
                u2 = self.tp2.real * (1 - eps)
                norbit = 1
            elif E < self.vmax and ust > self.tp3.real:
                u1 = 0.5
                u2 = self.tp3.real * (1 + eps)
                norbit = 0.5
            elif E > self.vmax:
                u1 = ust
                u2 = 0.5
                norbit = 0.5
        else:
            if E >= 0:
                u1 = ust
                u2 = 0.5
                norbit = 0.5
            else:
                u1 = self.tp1.real * (1 + eps)
                u2 = 0.5
                norbit = 0.5

        def theta(w):
            theta = (self.l / (2 ** (1 / 2))) * ((E - v(w, self.l)) ** (-1 / 2))
            return theta

        delphi, erro = quad(theta, u1, u2)

        n = 1000
        uc = np.arange(u1, u2, (u2 - u1) / n)
        ud = np.arange(u2, u1, (u1 - u2) / n)

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
            utotal = np.concatenate([uc, ud] * (norbit))

        accphi = [0] * (len(utotal))

        if norbit == 0.5:
            accphi = phi1
            x = [0] * (len(uc))
            y = [0] * (len(uc))
            for i in range(len(uc)):
                x[i] = (math.cos(accphi[i])) / utotal[i] * correction
                y[i] = (math.sin(accphi[i])) / utotal[i] * correction
        else:
            for i in range(norbit):
                for j in range(n):
                    accphi[j + (2 * i * n)] = 2 * i * delphi + phi1[j]
                    accphi[j + ((2 * i + 1) * n)] = ((2 * i) + 1) * delphi + phi2[j]
            x = [0] * (2 * norbit * n)
            y = [0] * (2 * norbit * n)
            for i in range(2 * norbit * n):
                x[i] = ((math.cos(accphi[i])) / utotal[i]) * correction
                y[i] = ((math.sin(accphi[i])) / utotal[i]) * correction

        #PLOT DA ÓRBITA

        fig2 = plt.figure()
        plt.plot(x, y, color="gold")
        plt.xlabel("x [km]")
        plt.ylabel("y [km]")
        circle = Circle((0, 0), 2 * 1.477, color='dimgrey')
        plt.gca().add_patch(circle)
        plt.gca().set_aspect('equal')
        plt.axis(
            [(-1 / u1 - 1) * 1.477, (1 / u1 + 1) * 1.477, (-1 / u1 - 1) * 1.477, (1 / u1 + 1) * 1.477])
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

    def calcularPlotPotencialVisaoTres(self):
        self.vlim = 0
        self.rmax = 80
        self.r = np.arange(2, self.rmax, self.rmax / 300)
        self.u = 1 / self.r

        #PLOT POTENCIAL

        fig1 = plt.figure()
        plt.subplot(1, 1, 1)
        plt.plot(1.477 * self.r, v(self.u, self.l), color="white")
        plt.axhline(0, linewidth=0.3, color='white')
        plt.xlabel("r [km]")
        plt.axis([0, 1.477 * self.rmax, -0.5, self.vlim + 0.1])
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

    def calcularPlotPotencialVisaoQuatro(self, energia, numeroOrbitas):
        print('energia=' + energia)
        print('numeroOrbitas=' + numeroOrbitas)

        E = eval(energia)
        norbit = eval(numeroOrbitas)

        if norbit < 1 or norbit > 20:
            print('Escala deve estar entre 1 e 20')
            return

        if E == 0:
            E = E + 1e-10
        
        self.coef = [- (self.l ** 2), (self.l ** 2) / 2, -1, -E]
        self.roots = np.roots(self.coef)
        self.tp1 = self.roots[2]
        self.tp2 = self.roots[1]
        self.tp3 = self.roots[0]

        eps = 0.00000001
        rst = 10 * self.l / (E + 0.5)
        ust = 1 / rst
        correction = 1.477

        if E >= 0:
            u1 = ust
            u2 = 0.5
            norbit = 0.5
        else:
            u1 = self.tp1.real * (1 + eps)
            u2 = 0.5
            norbit = 0.5

        def theta(w):
            theta = (self.l / (2 ** (1 / 2))) * ((E - v(w, self.l)) ** (-1 / 2))
            return theta

        delphi, erro = quad(theta, u1, u2)

        n = 1000
        uc = np.arange(u1, u2, (u2 - u1) / n)
        ud = np.arange(u2, u1, (u1 - u2) / n)

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
            utotal = np.concatenate([uc, ud] * (norbit))

        accphi = [0] * (len(utotal))

        if norbit == 0.5:
            accphi = phi1
            x = [0] * (len(uc))
            y = [0] * (len(uc))
            for i in range(len(uc)):
                x[i] = (math.cos(accphi[i])) / utotal[i] * correction
                y[i] = (math.sin(accphi[i])) / utotal[i] * correction
        else:
            for i in range(norbit):
                for j in range(n):
                    accphi[j + (2 * i * n)] = 2 * i * delphi + phi1[j]
                    accphi[j + ((2 * i + 1) * n)] = ((2 * i) + 1) * delphi + phi2[j]
            x = [0] * (2 * norbit * n)
            y = [0] * (2 * norbit * n)
            for i in range(2 * norbit * n):
                x[i] = ((math.cos(accphi[i])) / utotal[i]) * correction
                y[i] = ((math.sin(accphi[i])) / utotal[i]) * correction

        #PLOT ORBITA

        fig2 = plt.figure()
        plt.plot(x, y, color="gold")
        plt.xlabel("x [km]")
        plt.ylabel("y [km]")
        circle = Circle((0, 0), 2 * 1.477, color='dimgrey')
        plt.gca().add_patch(circle)
        plt.gca().set_aspect('equal')
        plt.axis([(-1 / u1 - 1) * 1.477, (1 / u1 + 1) * 1.477, (-1 / u1 - 1) * 1.477, (1 / u1 + 1) * 1.477])
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
