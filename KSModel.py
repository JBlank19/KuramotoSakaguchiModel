import os
import sys
import time
import matplotlib.pyplot as plt
import PyQt5
from threading import Thread

from matplotlib.animation import FuncAnimation
from matplotlib.widgets import TextBox
from pylab import *

from Distribuciones import *
from KuramotoSakaguchi import *

plt.rcParams['toolbar'] = 'None'


def initialize(K = 1, alpha = 1, T = 1, N = 100, num_en_pantalla = 10, mu = 0, sigma = 1):
    # Inicializa la distribución de frecuencias intrinsecas, la simulación (objeto), la lista de angulos iniciales
    # La lista que contiene el valor absoluto del parametro de orden, los indices aleatorios para el n en pantalla
    # la lista recortada con dichos indices y los colores que tendran los dots acorde a su omega natural
    global lista_orden, colors, muestra, solucion_activa, tiempo_activo, simulacion, random_indexes, omega_intrinseco

    dt = 2*np.pi / (sigma * 60)

    distribucion = Distribuciones(mu, sigma, N)
    omega_intrinseco = distribucion.Gaussiana()

    simulacion = KuramotoSakaguchi(K, alpha, dt, T, N, omega_intrinseco)
    theta_init = simulacion.theta_init()
    solucion_activa, tiempo_activo = simulacion.resolver_edo(theta_init, 0, T)
    # Filas son osciladores, columnas son tiempo
    lista_orden = np.array([])
    lista_orden = np.append(lista_orden, abs(simulacion.parametro_orden(theta_init)))

    random_indexes = np.random.randint(0, len(solucion_activa[:, 0]), num_en_pantalla)

    muestra=[]
    for z in range(num_en_pantalla):
        muestra.append(solucion_activa[random_indexes[z]][0])

    colors = omega_intrinseco[random_indexes]

    return K, alpha, T, dt, N, num_en_pantalla, mu, sigma, muestra, simulacion


def initialize_graphics():
    # Inicializa los ejes de matplotlib asociados a las gráficas
    global x_central, flecha_vertical, flecha_horizontal, cuarto_circulo_angulo, texto_phi, texto_radio
    fig = plt.figure()

    ax1 = plt.axes([0.03, 0.12, 0.32, 0.62]) # Eje del circulo con dots
    ax1.axis("equal")
    ax1.grid(False)
    ax1.axis('off')
    ax1.axes.set_title(r"Current system status", fontweight="bold", fontsize=13)
    ax2 = plt.axes([0.4, 0.1, 0.3, 0.3]) # Eje del gráfico de orden con el tiempo
    ax2.axes.set_title(r"Order parameter", fontweight="bold", fontsize=13)
    ax2.axes.set_ylabel(r"$|r e^{\phi i}|$", fontsize=13)

    ax3 = plt.axes([0.4, 0.52, 0.3, 0.3])  # Eje de la distribución normal

    # titulo superior
    axtitulo = plt.axes([0.25, 0.92, 0.5, 0.1])
    axtitulo.axis("off")
    axtitulo.text(0, 0, r"Kuramoto-Sakaguchi Model: $\dot{\theta}_i = \omega_i + \frac{K}{N} "
                        r"\sum_{j = 1}^N \sin{(\theta_j - \theta_i + \alpha)}, \ i = 1,...,N$", fontweight="bold", fontsize=15)

    rads = np.arange(0, (2 * np.pi), 0.01)
    ax1.plot(np.cos(rads), np.sin(rads), 'k', '-', zorder=0, lw=1)

    dots = ax1.scatter(np.cos(muestra), np.sin(muestra), c=colors, edgecolors="black", cmap="rainbow", s=150, zorder=1)
    dots.set_clim([min(omega_intrinseco), max(omega_intrinseco)])

    # Gráfico del circulo requiere ciertos objetos estáticos a pintar (ejes de r y phi, etc.)
    x_central = ax1.scatter(0, 0, c="black", marker="x")
    flecha_vertical = ax1.arrow(-1, -1, 0, 0.15, width=.01, facecolor='black', edgecolor='none')
    flecha_horizontal = ax1.arrow(-1, -1, 0.15, 0, width=.01, facecolor='black', edgecolor='none')

    x_circulo = np.linspace(-1, -0.9)
    y_circulo = np.sqrt(0.1**2 - (x_circulo + 1)**2) - 1
    cuarto_circulo_angulo, = ax1.plot(x_circulo, y_circulo, color="black")
    texto_phi = ax1.text(-0.92, -0.92, r"$\phi$", fontsize=13)
    texto_radio = ax1.text(-1.06, -0.84, r"r", fontsize=13)

    # Pintamos la linea de 1/sqrt(N) que será el valor de r para k = 0
    ax2.plot([0.01, 1.99], [1 / np.sqrt(N), 1 / np.sqrt(N)], lw=2, color="red", linestyle="dashed")
    ax2.text(1.8, 1/np.sqrt(N) + 0.03, r"$1 / \sqrt{N}$", color="red", fontsize=13)

    return fig, ax1, ax2, ax3, dots


def initialize_distribution(es_gaussiana = True):

    if es_gaussiana:
        global limite_x

        lista_opciones = [np.sqrt(N), 2]
        count, bins, paquetes = ax3.hist(omega_intrinseco, int(max(lista_opciones)), ec="black", density=True)
        limite_x = round(max(abs(bins)), 0)

        # Poner los bins de colores acorde a colormap
        cm = plt.cm.get_cmap('rainbow')
        bin_centers = 0.5 * (bins[:-1] + bins[1:])
        col = bin_centers - min(bin_centers)
        col /= max(col)

        for c, p in zip(col, paquetes):
            plt.setp(p, 'facecolor', cm(c))

        lista_x_gauss = np.linspace(-limite_x - 1, limite_x + 1, 100)
        ax3.plot(lista_x_gauss, 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(- (lista_x_gauss - mu) ** 2 / (2 * sigma ** 2)),
                  linewidth=2,
                  color='black')

        ax3.axes.set_xlim(-limite_x - 1, limite_x + 1)
        ax3.axes.set_title(r"Gaussian distribution of natural frequencies", fontweight="bold", fontsize=13)
        ax3.axes.set_ylabel(r"N (normalized)", fontsize=13)
        ax3.axes.set_xlabel(r"Natural frequencies [$\omega$]", fontsize=13)


        # plt.show()

def initialize_boxes():
    axtiempo = fig.add_axes([0.48, 0.05, 0.5, 0.05])
    axtiempo.axis("off")
    texto_time = axtiempo.text(0, 0.5, "", ha="left", va="top")

    # adjust the main plot to make room for the sliders
    # plt.subplots_adjust(right=0.7)

    width = 0.07
    height = 0.03

    # Axes for sliders / boxes / etc.
    ax_num_en_pantalla = plt.axes([0.87, 0.7, width, height])
    num_en_pantalla_box = TextBox(ax_num_en_pantalla, r'Number of osc. on screen  ', initial=num_en_pantalla, color="lavenderblush", hovercolor='thistle')
    ax_K = plt.axes([0.87, 0.6, width, height])
    K_box = TextBox(ax_K, r'Coupling constant [$K$]  ', initial=K, color="lavenderblush", hovercolor='thistle')

    ax_T = plt.axes([0.87, 0.4, width, height])
    T_box = TextBox(ax_T, r'Iteration period  ', initial=T, color="lightcyan", hovercolor='paleturquoise')
    ax_N = plt.axes([0.87, 0.5, width, height])
    N_box = TextBox(ax_N, r'Total number of oscillators [$N$]  ', initial=N, color="lightcyan", hovercolor='paleturquoise')
    ax_sigma = plt.axes([0.87, 0.3, width, height])
    sigma_box = TextBox(ax_sigma, r'Deviation [$\sigma$]  ', initial=sigma, color="lightcyan", hovercolor='paleturquoise')
    ax_alpha = plt.axes([0.87, 0.2, width, height])
    alpha_box = TextBox(ax_alpha, r'Alpha [$\alpha$]  ', initial=alpha, color="lightcyan", hovercolor='paleturquoise')

    return axtiempo, texto_time, ax_num_en_pantalla, num_en_pantalla_box, ax_K, K_box, ax_T, T_box, ax_N, N_box,\
           ax_sigma, sigma_box, ax_alpha, alpha_box

def initialize_buttons():
    width = 0.05
    height = 0.03

    # Axes botón reset
    ax_reset = plt.axes([0.87, 0.05, width, height])
    breset = Button(ax_reset, r"Reset", color="lightcyan", hovercolor='paleturquoise')
    breset.on_clicked(reset)

    # Axes botón cerrar
    ax_cerrar = plt.axes([0.81, 0.05, width, height])
    bcerrar = Button(ax_cerrar, r"Close", color="lightcyan", hovercolor="paleturquoise")
    bcerrar.on_clicked(cerrar)

    # Axes botón flecha interior
    ax_flecha_orden = plt.axes([0.87, 0.1, width, height])
    bflecha_orden = Button(ax_flecha_orden, r"$\checkmark$", color="lavenderblush", hovercolor='thistle')
    bflecha_orden.on_clicked(flecha_orden_show)
    ax_flecha_orden.text(-1, 0.2, r"Worm:")

    return ax_reset, breset, ax_flecha_orden, bflecha_orden, ax_cerrar, bcerrar


def initialize_slider():
    ax_slomo = plt.axes([0.87, 0.8, 0.07, 0.05])
    slomo_slider = Slider(
        ax=ax_slomo,
        label=r'Slow time  ',
        valmin=0,
        valmax=100,
        valstep=10,
        valinit=1,
        color="thistle",
        initcolor="none",
        track_color="lavenderblush",
        handle_style={'edgecolor':"black"}
    )
    return ax_slomo, slomo_slider


def reset(event):
    global anim
    global K, new_alpha, new_T, new_N, num_en_pantalla, mu, sigma

    anim.pause()
    plt.close()

    # Escribe en un archivo los parametros para reiniciar el script desde 0 con ellos.
    new_vars = [K, new_alpha, new_T, new_N, num_en_pantalla, mu, sigma]

    # open file in write mode
    with open(r'new_vars.txt', 'w') as fp:
        for item in new_vars:
            # write each item on a new line
            fp.write("%s\n" % item)
        print('Done writing variables to savefile.')

    print("Restarting!")
    os.execl(sys.executable, sys.executable, *sys.argv)


def cerrar(event):
    exit()

def flecha_orden_show(event):
    global show_flecha
    if show_flecha:
        show_flecha = False
        bflecha_orden.label.set_text(r"$\times$")
    else:
        show_flecha = True
        bflecha_orden.label.set_text(r"$\checkmark$")


def submit_num_en_pantalla(texto_num_en_pantalla):
    global num_en_pantalla, random_indexes, dots, muestra, colors
    num_en_pantalla = int(texto_num_en_pantalla)
    random_indexes = np.random.randint(0, len(solucion_activa[:, 0]), num_en_pantalla)

    # Dibujamos
    muestra = []
    for z in range(num_en_pantalla):
        muestra.append(solucion_activa[random_indexes[z]][-1])

    colors = omega_intrinseco[random_indexes]

    dots = ax1.scatter(np.cos(muestra), np.sin(muestra), c=colors, edgecolors="black", cmap="rainbow", s=150, zorder=1)


def submit_K(texto_K):
    global simulacion, K
    simulacion.update_K(float(texto_K))
    K = texto_K

def submit_alpha(texto_alpha):
    global new_alpha
    new_alpha = float(texto_alpha)

def submit_T(texto_T):
    global new_T
    new_T = float(texto_T)


def submit_N(texto_N):
    global new_N
    new_N = int(texto_N)


def submit_sigma(texto_sigma):
    global sigma
    sigma = float(texto_sigma)


# Inicializar variables desde el fichero.
variables_fichero = []

# Se abre el fichero creado
with open(r'new_vars.txt', 'r') as fp:
    for line in fp:
        x = line[:-1]
        variables_fichero.append(x)

# Se convierten a float o int las variables
for i in range(len(variables_fichero)):
    if i == 3 or i == 4 or i == 5:
        variables_fichero[i] = int(variables_fichero[i])
    else:
        variables_fichero[i] = float(variables_fichero[i])

# Se utilizan las variables preguardadas en el fichero para inicializar. La unica que se puede cambiar en directo
# posteriormente será K
K, alpha, T, dt, N, num_en_pantalla, mu, sigma, muestra, simulacion = initialize(*variables_fichero)
new_T, new_N, new_sigma, new_alpha = T, N, sigma, alpha
show_flecha = True

fig, ax1, ax2, ax3, dots = initialize_graphics()
initialize_distribution()

axtiempo, texto_time, ax_num_en_pantalla, num_en_pantalla_box, ax_K, K_box, ax_T, T_box, ax_N, N_box, \
ax_sigma, sigma_box, ax_alpha, alpha_box = initialize_boxes()
ax_reset, breset, ax_flecha_orden, bflecha_orden, ax_cerrar, bcerrar = initialize_buttons()
ax_slomo, slomo_slider = initialize_slider()

ax_credit = plt.axes([0, 0, 0.1, 0.1])
ax_credit.axis("off")
ax_credit.text(0.1, 0.1, r"Made by J. Blanco and E. del Campo", fontsize=10)

line_orden, = ax2.plot(tiempo_activo[0], lista_orden[0], lw=2, color="black")
flecha_orden, = ax1.plot(lista_orden[0], lista_orden[0], lw=2, solid_capstyle='round', color=str(line_orden.get_color()))

# register the update function with each slider
num_en_pantalla_box.on_submit(submit_num_en_pantalla)
K_box.on_submit(submit_K)
T_box.on_submit(submit_T)
N_box.on_submit(submit_N)
sigma_box.on_submit(submit_sigma)
alpha_box.on_submit(submit_alpha)

i_stamp = 0
results = {}
calculation = None
tiempo_total = tiempo_activo


def calculate(prev_sol, start_time, end_time, result):
    # Aqui result es un objeto mutable (ya que los Threads no pueden devolver)
    result['sols'], result['time'] = simulacion.resolver_edo(prev_sol, start_time, end_time)
    #print("Thread de calculo ha finalizado")


def animate(i):
    global results, solucion_activa, tiempo_activo, lista_orden, tiempo_total, T
    global i_stamp, calculation
    global simulacion, anim
    global colors

    try:
        # The first frame of cycle
        if i - i_stamp == 1 and calculation is None:
            #print("Iniciamos thread de calculo")
            results = {}
            calculation = Thread(target=calculate,
                                 args=(solucion_activa[:, -1], tiempo_activo[-1], tiempo_activo[-1] + T, results))
            calculation.start()

        # Last frame of cycle
        elif len(solucion_activa[0, :]) - 1 == i - i_stamp:
            #print("Uniendo threads")
            calculation.join()
            solucion_activa = results['sols']
            tiempo_activo = results['time']
            tiempo_total = np.append(tiempo_total, tiempo_activo)
            calculation = None
            i_stamp = i

        # Dibujamos
        muestra = []
        for z in range(num_en_pantalla):
            muestra.append(solucion_activa[random_indexes[z]][i - i_stamp])

    except Exception as e:
        print(e)
        exit(1)

    dots.set_offsets(np.c_[np.cos(muestra), np.sin(muestra)])

    orden_actual = simulacion.parametro_orden(solucion_activa[:, i - i_stamp])

    lista_orden = np.append(lista_orden, orden_actual)
    abs_lista_orden = abs(lista_orden)
    ang_lista_orden = np.arctan2(lista_orden.imag, lista_orden.real)
    line_orden.set_data(tiempo_total[:i], abs_lista_orden[:i])

    if show_flecha:
        flecha_orden.set_visible(True)
        x_central.set_visible(True)
        flecha_vertical.set_visible(True)
        flecha_horizontal.set_visible(True)
        cuarto_circulo_angulo.set_visible(True)
        texto_phi.set_visible(True)
        texto_radio.set_visible(True)

        flecha_orden.set_data(abs_lista_orden[-20:i] * np.cos(ang_lista_orden[-20:i]),
                              abs_lista_orden[-20:i] * np.sin(ang_lista_orden[-20:i]))
    elif not show_flecha:
        flecha_orden.set_visible(False)
        x_central.set_visible(False)
        flecha_vertical.set_visible(False)
        flecha_horizontal.set_visible(False)
        cuarto_circulo_angulo.set_visible(False)
        texto_phi.set_visible(False)
        texto_radio.set_visible(False)

    ax2.axis([0, tiempo_total[-1] + tiempo_total[-1] / 5, 0, 1])
    ax2.set_xticks([0, 2])
    ax2.set_xticklabels([0, r"$t_{max}$"])

    texto_time.set_text("Elapsed time: " + str(round(tiempo_activo[i - i_stamp], 1)) +
                        " of " + str(int(tiempo_activo[-1])) + " steps")

    time.sleep(slomo_slider.val/1000)

    return dots, texto_time, line_orden, flecha_orden, x_central, flecha_vertical, flecha_horizontal, \
           cuarto_circulo_angulo, texto_phi, texto_radio


anim = FuncAnimation(fig, animate, interval=1, blit=True, repeat=False)

figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
figManager.set_window_title('Kuramoto-Sakaguchi Model')
plt.show()

__all__ = ["os", "sys", "time", "Distribuciones"]
