import numpy as np


def solver_rk4(f, a, b, r0, N = 1000, args=()):
    # Calcularemos la solucion dr/dt=f(t,r) con valor inicial r(a)=r0, y con paso h=(b-a)/N, usando RK4.
    # Notese que r0 tiene que ser un vector, esto es, un np.array, en cuyo caso r(t) sera una funcion vectorial,
    # de la que nos devolvera el metodo: r(0), r(h), r(2h),...r(Nh=b).
    # Por lo tanto, lo que nos devuelve el metodo es una matriz "T" con N+1 filas y
    # numero de columnas igual a la dimension de r0.


    # Fijamos el tamaño de sol (Columnas: osciladores, filas: tiempo)
    sol = np.zeros((N+1, r0.size))

    # Se escribe la primera fila de la solución como el valor inicial dado a la función
    sol[0, :] = r0

    # Calculamos el valor de h
    h = (b-a)/N

    # Hacemos el bucle de RK4
    for k in range(1, N + 1):
        k1 = h * f((k - 1) * h, sol[k - 1, :], *args)
        k2 = h * f((k - 1 / 2) * h, sol[k - 1, :] + k1 / 2, *args)
        k3 = h * f((k - 1 / 2) * h, sol[k - 1, :] + k2 / 2, *args)
        k4 = h * f((k) * h, sol[k - 1, :] + k3, *args)
        sol[k, :] = sol[k - 1, :] + (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
    return sol.T

