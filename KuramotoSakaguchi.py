import RK4
from Distribuciones import *
# from scipy.integrate import odeint
from RK4 import *
import cmath


class KuramotoSakaguchi:
    def __init__(self, K, alpha, dt, T, N, lista_omega):
        self.K = K
        self.alpha = alpha
        self.dt = dt
        self.T = T
        self.N = N
        self.lista_omega = lista_omega
        self.r = 0
        self.phi = 0

    def theta_init(self):
        return 2 * np.pi * np.random.random_sample(self.N)

    def generador_edo(self, lista_t, lista_theta, acoplamiento):
        """
        Recibe una lista de thetas y de omegas intrinsecas y genera una lista de dthetadt de acuerdo a Kuramoto.
        """
        comp = self.parametro_orden(lista_theta)
        self.r, self.phi = abs(comp), cmath.phase(comp)

        # angulo_i, angulo_j = np.meshgrid(lista_theta, lista_theta)
        # intercambio = np.sin(angulo_j - angulo_i)

        dthetadt = self.lista_omega + acoplamiento * self.r * np.sin(self.phi - lista_theta + self.alpha)

        return dthetadt

    def resolver_edo(self, lista_theta, a, b):

        lista_t = np.linspace(a, b, int(((b-a) / (self.dt))))
        # sol = odeint(self.generador_edo, lista_theta, lista_t, args=(acoplamiento,))
        sol = RK4.solver_rk4(self.generador_edo, a, b, lista_theta, int(self.T / self.dt), args=(self.K,))
        return sol, lista_t

    def parametro_orden(self, lista_theta):
        # Devuelve r y phi del parametro de orden
        complejo = (1 / self.N) * sum([(np.e ** (theta * 1j)) for theta in lista_theta])

        return complejo

    def update_K(self, new_K):
        self.K = new_K

    def update_alpha(self, new_alpha):
        self.alpha = new_alpha

    def update_T(self, new_T):
        self.T = new_T

    def update_dt(self, new_dt):
        self.dt = new_dt

    def update_N(self, new_N):
        self.N = new_N
