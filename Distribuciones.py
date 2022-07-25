import numpy as np
import matplotlib.pyplot as plt


class Distribuciones:
    def __init__(self, mu, sigma, N):
        self.mu = mu
        self.sigma = sigma
        self.N = N

    def Gaussiana(self):
        """
        Devuelve una lista de omegas intrinsecas según una gaussiana
        Input: media, desv std, numero
        Output: lista con valores omega
        """
        lista_gauss = np.random.normal(self.mu, self.sigma, self.N)

        return lista_gauss
