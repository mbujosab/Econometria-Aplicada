import numpy as np
import pandas as pd
import matplotlib as mpl
# definimos parámetros para mejorar los gráficos
mpl.rc('text', usetex=False)
import matplotlib.pyplot as plt

import statsmodels.api as sm
from statsmodels.tsa.arima_process import ArmaProcess
from scipy import signal
from scipy.signal import freqz
from numpy.polynomial import Polynomial
import math

def theoretical_pacf_from_acf(acf_vals):
    """
    Calcula la PACF teórica (para lags >= 0) de forma paramétrica usando el algoritmo de Durbin-Levinson.

    Parámetros:
      acf_vals (list):
          una lista tal que acf_vals[0] corresponde al retardo 0.
    
    Retorna:
      list: el coeficiente en la posición k corresponde a la PACF en el retardo k.
    """
    nlags = len(acf_vals)-1
    # Matriz para almacenar los coeficientes (phi)
    phi = np.zeros((nlags+1, nlags+1))
    # v guarda la varianza de error en cada paso
    v = np.zeros(nlags+1)
    pacf = np.zeros(nlags+1)
    v[0] = acf_vals[0]
    pacf[0] = 1.0  # Por convención, PACF(0)=1.
    
    for k in range(1, nlags+1):
        if k == 1:
            phi[1,1] = acf_vals[1] / acf_vals[0]
            pacf[1] = phi[1,1]
            v[1] = v[0]*(1 - phi[1,1]**2)
        else:
            sum_val = 0.0
            # Sumar phi[k-1,j]*acf_vals[k-j] para j=1,...,k-1.
            for j in range(1, k):
                sum_val += phi[k-1, j] * acf_vals[k - j]
            phi[k, k] = (acf_vals[k] - sum_val) / v[k-1]
            pacf[k] = phi[k, k]
            for j in range(1, k):
                phi[k, j] = phi[k-1, j] - phi[k, k] * phi[k-1, k - j]
            v[k] = v[k-1] * (1 - phi[k, k]**2)
    return pacf  # Incluye PACF en lag 0 con valor 1.

def plot_arma_parametric_diagnostics(ar_params=[1,], ma_params=[1,], sigma2=1, lags=20):
    """
    Genera y retorna una figura con tres subgráficas para un modelo ARMA:
      - ACF teórica (omitiendo retardo 0)
      - PACF teórica paramétrica (calculada vía el algoritmo de Durbin-Levinson, omitiendo retardo 0)
      - Densidad espectral teórica
    
    Parámetros:
      ar_params (list):
          Coeficientes del polinomio AR (debe incluir el 1 inicial).
      ma_params (list):
          Coeficientes del polinomio MA (debe incluir el 1 inicial).
      sigma2 (float):
          Varianza del ruido Uₜ.
      lags (int, opcional): 
          Número de retardos considerados (por defecto 20).
    
    Retorna:
      fig : objeto matplotlib.figure.Figure
    """
    
    phi   = Polynomial(ar_params, symbol='B')
    theta = Polynomial(ma_params, symbol='B')
    # Crear el proceso ARMA
    arma_process = ArmaProcess(ar=phi.coef, ma=theta.coef)
    
    # 1. Calcular la ACF teórica
    acf_theo = arma_process.acf(lags=lags)
    lags_theo = np.arange(len(acf_theo))
    # Omitir el retardo 0 para gráficos de ACF y PACF:
    acf_plot = acf_theo[1:]
    lags_plot = lags_theo[1:]
    
    # 2. Calcular la PACF teórica de forma paramétrica a partir de la ACF
    pacf_theo = theoretical_pacf_from_acf(acf_theo)
    pacf_plot = pacf_theo[1:]
    lags_pacf = np.arange(len(pacf_theo))[1:]
    
    # 3. Calcular la densidad espectral teórica usando freqz
    w, h = freqz(theta.coef, phi.coef, worN=1024)
    freq = w / (2 * np.pi)  # Conversión de radianes/muestra a ciclos/muestra
    spectrum = (sigma2 / (2 * np.pi)) * np.abs(h)**2
    
    # 4. Graficar los tres subgráficos
    fig, axs = plt.subplots(1, 3, figsize=(18, 3))
    
    # ACF teórica
    axs[0].stem(lags_plot, acf_plot, basefmt=" ")
    axs[0].axhline(0, color='black', lw=0.5, linestyle='--')
    axs[0].set_xlabel('Retardo')
    axs[0].set_ylabel('ACF')
    
    # PACF teórica paramétrica
    axs[1].stem(lags_pacf, pacf_plot, basefmt=" ")
    axs[1].axhline(0, color='black', lw=0.5, linestyle='--')
    axs[1].set_xlabel('Retardo')
    axs[1].set_ylabel('PACF')
    
    # Densidad espectral teórica
    axs[2].plot(freq, spectrum)
    axs[2].axhline(0, color='black', lw=0.5, linestyle='--')
    axs[2].set_xlabel('Frecuencia')
    axs[2].set_ylabel('Densidad Espectral')
    
    plt.tight_layout()
    plt.close(fig)
    return fig

def plot_arma_analysis(ar_params=[1,], ma_params=[1,], sigma2=1, lags=20, n=400, seed=0):
    """
    Simula una serie temporal utilizando un modelo ARMA y visualiza la serie generada junto con sus funciones de autocorrelación (ACF), autocorrelación parcial (PACF) y el periodograma estimado.
    
    Parameters
    ----------
    ar_params : list, optional
        Coeficientes AR del modelo (por defecto [1]).
    ma_params : list, optional
        Coeficientes MA del modelo (por defecto [1]).
    sigma2 : float, optional
        Varianza del ruido blanco (por defecto 1).
    lags : int, optional
        Número de lags a incluir en ACF y PACF (por defecto 20).
    n : int, optional
        Número de observaciones a simular (por defecto 400).
    seed : int, optional
        Semilla para la generación aleatoria (por defecto 0).
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        La figura generada que incluye la serie temporal y los gráficos ACF, PACF y el periodograma.
    """

    lags=math.floor(n/2) if n/2 < lags else lags
    
    phi   = Polynomial(ar_params, symbol='B')
    theta = Polynomial(ma_params, symbol='B')
    arma_process = sm.tsa.ArmaProcess(phi.coef, theta.coef)
    
    # Simulación de la serie temporal
    if seed:
        np.random.seed(seed)
    data = arma_process.generate_sample(nsample=n)
       
    # Crear la figura
    fig = plt.figure(figsize=(18, 7))
    
    # Subgráfico de la serie temporal
    ax1 = fig.add_subplot(2, 1, 1)
    ax1.plot(data, color='blue')
    ax1.set_xlabel('t')
    
    # Crear una fila de 3 subgráficos para ACF, PACF y periodograma
    ax2 = fig.add_subplot(2, 3, 4)  # Fila 2, Columna 1
    sm.graphics.tsa.plot_acf(data, lags=lags, ax=ax2, title='ACF Estimada', zero=False, auto_ylims=True)
    ax2.set_xlabel('Retardo')

    ax3 = fig.add_subplot(2, 3, 5)  # Fila 2, Columna 2
    sm.graphics.tsa.plot_pacf(data, lags=lags, ax=ax3, title='PACF Estimada', zero=False, auto_ylims=True)
    ax3.set_xlabel('Retardo')
    
    ax4 = fig.add_subplot(2, 3, 6)  # Fila 2, Columna 3
    f, Pxx = signal.welch(data, fs=1, nperseg=256)
    ax4.plot(f, Pxx, color='purple')
    ax4.set_title('Periodograma Estimado')
    ax4.set_xlabel('Frecuencia')
    #ax4.set_ylabel('Densidad Espectral')
    ax4.set_xlim(0, 0.5)  # Limitar el eje x a 0.5

    plt.tight_layout()
    plt.close(fig)    
    return fig

def polynomial_roots_table(coef_polinomio):
    """
    Función que calcula las raíces de un polinomio y presenta sus características en forma de tabla.

    Args:
        coef_polinomio (list): Lista de coeficientes del polinomio, donde el índice representa el grado del término.

    Returns:
        pd.DataFrame: DataFrame que contiene las raíces del polinomio junto a su parte real, parte imaginaria, 
                      módulo, frecuencia, período y multiplicidad.
    
    Descripción:
        Esta función utiliza la librería SymPy para obtener las raíces de un polinomio definido por sus coeficientes. 
        Luego, calcula la parte real e imaginaria de cada raíz, así como su módulo, frecuencia correspondiente 
        a la fase de la raíz, período asociado a la frecuencia y la multiplicidad de cada raíz. Finalmente, 
        se presenta la información en un DataFrame de pandas.
    """
    polinomio = Polynomial(coef_polinomio, symbol='B')
    roots = polinomio.roots()
    
    # Contar la multiplicidad de las raíces
    unique_roots, counts = np.unique(roots, return_counts=True)

    # Crear listas para cada columna de la tabla
    real_parts = []
    imaginary_parts = []
    moduli = []
    frequencies = []
    periods = []
    multipliers = []

    for root, count in zip(unique_roots, counts):
        real = root.real
        imag = root.imag
        
        # Calcular módulo
        modulus = np.abs(root)
        moduli.append(modulus)

        # Calcular frecuencia y período
        frequency = np.abs(np.angle(root)) / (2*np.pi)
        period = 1/frequency if frequency else 'Inf'
        
        # Append results
        real_parts.append(real)
        imaginary_parts.append(imag)
        frequencies.append(frequency)
        periods.append(period)
        multipliers.append(count)  # Multiplicidad

    # Crear un DataFrame para mostrar la tabla
    table = pd.DataFrame({
        'P. real': real_parts,
        'P. imag.': imaginary_parts,
        'Módulo': moduli,
        'Frecuencia': frequencies,
        'Periodo': periods,
        'Multip.': multipliers
    })
    
    # Ordenar la tabla por la columna 'Frequency' en orden ascendente
    tabla = table.sort_values(by='Frecuencia').round(6)   
    return tabla

def Localizacion_parametros_plot(parametro1, parametro2, modelo='MA'):
    """
    Genera un gráfico que localiza el par de parámetros de un polinomio de grado 2 en relación con un triángulo,
    indicando las regiones donde dicho polinomio es invertible.

    Args:
        parametro1 (float): Valor del primer parámetro del polinomio.
        parametro2 (float): Valor del segundo parámetro del polinomio.
        modelo (str, opcional): Tipo de modelo a utilizar. Por defecto es 'MA'. Puede ser 'AR' para ajustar las 
                                representaciones de los parámetros.

    Returns:
        matplotlib.figure.Figure: Objeto Figure que contiene el gráfico generado.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_aspect('equal')

    # Set axis
    xmin, ymin = -2.2, -1.2
    xmax, ymax = -xmin, -ymin
    plt.axis([xmin, xmax, ymin, ymax])

    # Define the Greek letter based on modelo with subscripts
    simb_param_1 = r'$\phi_1$' if modelo == 'AR' else r'$\theta_1$'
    simb_param_2 = r'$\phi_2$' if modelo == 'AR' else r'$\theta_2$'
    
    # Set axis labels
    ax.set(xticks=[], yticks=[])
    ax.set_xlabel(simb_param_2, fontsize=16)
    ax.xaxis.set_label_position('top')
    ax.set_ylabel(simb_param_1, rotation=0, fontsize=16)
    ax.yaxis.set_label_position('right')
    
    # Draw (t1, t2) points
    p1 = np.linspace(-2, 2, 100)
    ax.plot(p1, -abs(p1) + 1, c='black', linestyle='--')
    ax.plot(p1, np.full_like(p1, -1), c='lightgray', linestyle='--')
    ax.plot(p1, -(p1**2 / 4), c='black')
    
    # Fill the region below the parabola and above the base of the triangle
    ax.fill_between(p1, -(p1**2 / 4), -1, where=-(p1**2 / 4) >= -1, color='lightgray', alpha=0.5)
    
    # Turn normal axes off
    for spine in ['left', 'bottom', 'top', 'right']:
        ax.spines[spine].set_visible(False)
        
    # Add arrows to represent axes
    axes_arrows = {'arrowstyle': '<|-|>', 'lw': 0.5, 'color': 'gray'}  # Hacerlas grises y más finas
    ax.annotate('', xy=(xmin, 0), xytext=(xmax, 0), arrowprops=axes_arrows)
    ax.annotate('', xy=(0, ymin), xytext=(0, ymax), arrowprops=axes_arrows)
    
    # Label complex region
    ax.text(0.09, -0.8, 'Raíces complejas', ha='center', fontsize=14, color='gray')
    

    # Add small marker to y-axis
    ax.axhline(y=1.005, xmin=0.495, xmax=0.505, c='black')
    ax.text(-0.12, -1.12, '-1', fontsize=10)
    ax.text(-0.12, 0.98, '1', fontsize=10)

    # Add small markers to x-axis
    ax.axvline(x=-2, ymin=0.495, ymax=0.505, c='black')
    ax.axvline(x=2, ymin=0.495, ymax=0.505, c='black')
    ax.text(-2.12, -0.12, '-2', fontsize=10)
    ax.text(1.88, -0.12, '2', fontsize=10)

    marker = ax.scatter(parametro1, parametro2, s=100, color='blue', marker='.')
    ax.legend([marker], [f'({simb_param_1}={parametro1}, {simb_param_2}={parametro2})'], loc='upper right', fontsize=12)
    plt.tight_layout()
    plt.close(fig)
    
    return fig

def RegionParametrica(modelo='generico'):   
    """
    Dibuja y anota un gráfico que representa regiones paramétricas dependiendo del modelo especificado.

    Args:
        modelo (str): Tipo de modelo a utilizar. Puede ser 'AR', 'MA' o 'generico'. 
                      'AR' representa un modelo autorregresivo, 'MA' un modelo de medias móviles y 
                      cualquier otro valor utilizará un modelo genérico.

    Returns:
        figura: Un objeto figura que representa el gráfico generado.
    
    Ejemplo:
        >>> fig = RegionParametrica(modelo='AR')
        >>> plt.show(fig)
    """  

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_aspect('equal')

    # Set axis
    xmin, ymin = -3, -2
    xmax, ymax = -xmin, -ymin
    plt.axis([xmin, xmax, ymin, ymax])

    # Define the Greek letter based on modelo with subscripts
    if modelo == 'AR':
        simb_param_1 = r'$\phi_1$'
        simb_param_2 = r'$\phi_2$'
    elif modelo == 'MA':
        simb_param_1 = r'$\theta_1$'
        simb_param_2 = r'$\theta_2$'
    else:
        simb_param_1 = r'$\tau_1$'
        simb_param_2 = r'$\tau_2$'

    desdolariza = lambda texto: texto.replace('$', '')

    # Set axis labels
    ax.set(xticks=[], yticks=[])
    ax.set_xlabel(simb_param_2, fontsize=16)
    ax.xaxis.set_label_position('top')
    ax.set_ylabel(simb_param_1, rotation=0, fontsize=16)
    ax.yaxis.set_label_position('right')
    
    # Draw (t1, t2) points
    p1 = np.linspace(-2, 2, 100)
    ax.plot(p1, -abs(p1) + 1, c='black', linestyle='--')
    ax.plot(p1, np.full_like(p1, -1), c='lightgray', linestyle='--')
    ax.plot(p1, -(p1**2 / 4), c='black')
    
    # Fill the region below the parabola and above the base of the triangle
    ax.fill_between(p1, -(p1**2 / 4), -1, where=-(p1**2 / 4) >= -1, color='lightgray', alpha=0.5)
    
    # Turn normal axes off
    for spine in ['left', 'bottom', 'top', 'right']:
        ax.spines[spine].set_visible(False)
        
    # Add arrows to represent axes
    axes_arrows = {'arrowstyle': '<|-|>', 'lw': 0.5, 'color': 'gray'}  # Hacerlas grises y más finas
    ax.annotate('', xy=(xmin, 0), xytext=(xmax, 0), arrowprops=axes_arrows)
    ax.annotate('', xy=(0, ymin), xytext=(0, ymax), arrowprops=axes_arrows)
    
    # Label complex region
    ax.text(0.09, -0.8, 'Raíces complejas', ha='center', fontsize=14, color='gray')
    

    # Annotate the plot with equations
    plot_arrowsl = {'arrowstyle': '-|>', 'connectionstyle': "arc3, rad=-0.2"}
    plot_arrowsr = {'arrowstyle': '-|>', 'connectionstyle': "arc3, rad=0.2"}
    ax.annotate('$' + desdolariza(simb_param_1) + '+' + desdolariza(simb_param_2) + '<1$', xy=(0.4, 0.25), xytext=(1.1, 0.4), arrowprops=plot_arrowsr, fontsize='12')
    ax.annotate('$' + desdolariza(simb_param_1) + '+' + desdolariza(simb_param_2) + '=1$', xy=(0.38, 0.6), xytext=(0.6, 0.8), arrowprops=plot_arrowsr, fontsize='12')
    
    ax.annotate('$' + desdolariza(simb_param_2) + '<1+' + desdolariza(simb_param_1)+ '$', xy=(-0.4, 0.25), xytext=(-1.7, 0.4), arrowprops=plot_arrowsl, fontsize='12')
    ax.annotate('$' + desdolariza(simb_param_2) + '=1+' + desdolariza(simb_param_1)+ '$', xy=(-0.38, 0.6), xytext=(-1.2, 0.8),arrowprops=plot_arrowsl, fontsize='12')

    ax.annotate('$' + desdolariza(simb_param_2) + '=-1$', xy=(1.5, -1), xytext=(1.8, -1.3), arrowprops=plot_arrowsl, fontsize='12')

    ax.annotate('$' + desdolariza(simb_param_1) + '^2 + 4' + desdolariza(simb_param_2) + '= 0$', xy=(1.0, -0.28), xytext=(1.8, -0.5), arrowprops=plot_arrowsr, fontsize='12')
    ax.annotate('$' + desdolariza(simb_param_1) + '^2 + 4' + desdolariza(simb_param_2) + '< 0$', xy=(-.9, -0.6), xytext=(-2.8, -0.5), arrowprops=plot_arrowsl, fontsize='12')
   
    # Label cases
    ax.text(.3, 0.25, '1', ha='center', fontsize=16, color='blue')
    ax.text(-.3, 0.25, '2', ha='center', fontsize=16, color='blue')
    ax.text(-.3, -0.35, '3', ha='center', fontsize=16, color='blue')
    ax.text(.3, -0.35, '4', ha='center', fontsize=16, color='blue')

    # Add small marker to y-axis
    ax.axhline(y=1.005, xmin=0.495, xmax=0.505, c='black')
    ax.text(-0.12, -1.12, '-1', fontsize=10)
    ax.text(-0.12, 0.98, '1', fontsize=10)

    # Add small markers to x-axis
    ax.axvline(x=-2, ymin=0.495, ymax=0.505, c='black')
    ax.axvline(x=2, ymin=0.495, ymax=0.505, c='black')
    ax.text(-2.12, -0.12, '-2', fontsize=10)
    ax.text(1.88, -0.12, '2', fontsize=10)
    
    plt.tight_layout()
    plt.close(fig)    
    
    return fig

def RangoMediana(datos, k=6, ax=None):
    # Inicializar listas para almacenar las medianas y rangos
    medians = []
    ranges = []
    
    n = len(datos)
    sub_sample_size = n // k
    for i in range(k):
        start = i * sub_sample_size
        end = start + sub_sample_size
        
        if i == k - 1:  # Para el último, incluir los datos sobrantes
            end = n
        
        sub_sample = datos[start:end]
        
        # Calcular mediana y rango
        median = np.median(sub_sample)
        range_value = np.max(sub_sample) - np.min(sub_sample)
        
        medians.append(median)
        ranges.append(range_value)

    # Crear el gráfico rango-mediana en el eje especificado
    if ax is not None:
        ax.scatter(medians, ranges) 
        ax.set_title('Gráfico Rango-Mediana')
        ax.set_xlabel('Mediana de Submuestra')
        ax.set_ylabel('Rango de Submuestra (Max - Min)')
        ax.grid()
        
# def RangoMediana(datos, k=6):
#     # Inicializar listas para almacenar las medianas y rangos
#     medians = []
#     ranges = []
#     
#     # Calcular el tamaño de las submuestras
#     n = len(datos)
#     sub_sample_size = n // k
#     for i in range(k):
#         start = i * sub_sample_size
#         end = start + sub_sample_size
#         
#         if i == k - 1:  # Para el último, incluir los datos sobrantes
#             end = n
#         
#         sub_sample = datos[start:end]
#         
#         # Calcular mediana y rango
#         median = np.median(sub_sample)
#         range_value = np.max(sub_sample) - np.min(sub_sample)
#         
#         medians.append(median)
#         ranges.append(range_value)
# 
#     # Crear el gráfico rango-mediana
#     fig, ax = plt.subplots(figsize=(5, 4))
#     ax.scatter(medians, ranges) # , color='r')
#     ax.set_title('Gráfico Rango-Mediana')
#     ax.set_xlabel('Mediana de Submuestra')
#     ax.set_ylabel('Rango de Submuestra (Max - Min)')
#     ax.grid()
#     plt.tight_layout()
#     plt.close(fig)    
#     
#     return fig
#
