from sympy import latex, Rational
from fractions import Fraction
import sympy
from sympy import Add
import matplotlib.pyplot as plt
from IPython.display import display, Math, display_png, Latex
import tempfile
from os.path import join           

def html(TeX):
    """ Plantilla HTML para insertar comandos LaTeX """
    return "<p style=\"text-align:center;\">$" + TeX + "$</p>"

def latex(a):
    """Método latex general"""
    try:
        return a.latex()
    except:
        return sympy.latex(a)

def pinta(data):
    """Muestra en Jupyter la representación latex de data"""
    display(Latex('$$'+latex(data)+'$$'))
    #display(Math(latex(data)))
    
class SerieConPrincipio:
    """
    Clase para representar series con un índice inicial arbitrario.

    Atributos:
    - coeficientes (list): Lista de coeficientes de la serie, sin ceros iniciales.
    - cogrado (int): Índice del primer coeficiente no nulo de la serie.
    - usa_float (bool): Indica si la serie opera con números flotantes.
    - final (bool): Indica si el final de la serie se representa con una elipsis.

    Métodos:
    - __repr__: Representa la serie como un par (cogrado, lista de coeficientes).
    - _repr_latex_: Representa la serie en formato LaTeX.
    - __add__: Suma dos SerieConPrincipio.
    - __mul__: Producto por un escalar o producto convolución entre SerieConPrincipio.
    - __rmul__: Producto por un escalar desde la izquierda.
    - inversa: Calcula la SerieConPrincipio inversa.
    """
    def __init__(self, coeficientes, cogrado=0, final=False, RepEscalera=False):
        primeros_no_nulos = next((i for i, c in enumerate(coeficientes) if c != 0), None)
        if primeros_no_nulos is None:
            self.coeficientes = []
            self.cogrado = 0
        else:
            self.coeficientes = coeficientes[primeros_no_nulos:]
            self.cogrado = cogrado + primeros_no_nulos
    
        self.final = final
        self.usa_float = any(isinstance(c, float) for c in self.coeficientes if not hasattr(c, 'is_number'))
        self.coeficientes = [
            float(c) if self.usa_float and isinstance(c, (int, float)) else Fraction(c) if isinstance(c, (int, float)) else c
            for c in self.coeficientes
        ]
        self.RepEscalera=RepEscalera
    
    def __repr__(self):
        """
        Representa la serie como un par (cogrado, lista de coeficientes).
    
        Retorna:
        - str: Representación de la serie como un par.
        """
        return f"({self.cogrado}, {self.coeficientes})"
    
    def latex(self):
        """
        Representa la serie en formato LaTeX.
    
        Retorna:
        - str: Representación en formato LaTeX.
        """
    
        if not self.coeficientes:
            return "$0$"
    
        terminos = []
        for i, coef in enumerate(self.coeficientes):
            indice = self.cogrado + i
            if coef != 0:
                coef_sympy = Rational(coef.numerator, coef.denominator) if isinstance(coef, Fraction) else coef
                coef_latex = r"\left(" + latex(coef_sympy) + r"\right)" if  isinstance(coef_sympy, Add) else latex(coef_sympy)
                if coef_sympy == -1 and indice != 0:
                    coef_latex = "-"
                if coef_sympy == 1 and indice != 0:
                    coef_latex = ""
                if indice == 0:
                    terminos.append(f"{coef_latex}")
                elif indice == 1:
                    terminos.append(f"{coef_latex}z")
                else:
                    terminos.append(f"{coef_latex}z^{{{indice}}}")  # Aseguramos que el exponente esté entre llaves
        cadena = r' + \cdots' if self.final else ''
        #return '$'+ " + ".join(terminos).replace("+ -", "- ") + cadena +'$'
        return " + ".join(terminos).replace("+ -", "- ") + cadena 
    
    def _repr_latex_(self):
        """
        Representa la serie en formato LaTeX.
    
        Retorna:
        - str: Representación en formato LaTeX para notebooks de Jupyter.
        """
        return '$$'+self.latex()+'$$'
    
    def _repr_html_(self):
        """ Construye la representación para el entorno jupyter notebook """
        if self.RepEscalera:
            return html(self.escalera())
        else:
            return html(self.latex())
    
    #def _repr_png_(self):
    #    """ Representación png para el entorno jupyter en Emacs """
    #    try:
    #        expr = '$'+self.latex()+'$'
    #        workdir = tempfile.mkdtemp()
    #        with open(join(workdir, 'borrame.png'), 'wb') as outputfile:
    #            sympy.preview(expr, viewer='BytesIO', outputbuffer=outputfile)
    #        return open(join(workdir, 'borrame.png'),'rb').read()
    #    except:
    #        return '$'+self.latex()+'$'
    
    def coeficientesEscalera(self):
        coeficientes = list(self.coeficientes)  # Creamos una copia de la lista de coeficientes
    
        if self.cogrado > 0:
            coeficientes = [0]*self.cogrado + coeficientes
            origen = 0        
        else: # len(self.coeficientes) + self.cogrado > 0:
            coeficientes = coeficientes + [0]*(1 - self.cogrado - len(self.coeficientes))
            origen = -self.cogrado
            
        return coeficientes, origen
    
    def escalera(self):
        coeficientes, origen = self.coeficientesEscalera()
        celdas = [r'\cdots 0'] + [r'{\color{blue}{\{' + str(c) + '\}}}' if i == origen else f"{str(c)}" for i,c in enumerate(coeficientes)] + [r'\cdots']    
        size = len(celdas)-1
        contenido = ' & '.join(celdas)
        latex_representation = "\\begin{array}{" + "c|"*size + "c}\n\\hline\n"
        latex_representation += f"{contenido} \\\\\hline"
        latex_representation += r"\end{array}"
        return '$' + latex_representation + '$'
    
    def __add__(self, otra):
        """
        Suma dos series.
    
        Parámetros:
        - otra (SerieConPrincipio): Otra serie.
    
        Retorna:
        - SerieConPrincipio: La serie resultante de la suma.
        """
        cogrado = min(self.cogrado, otra.cogrado)
        indice_final = max(
            self.cogrado + len(self.coeficientes),
            otra.cogrado + len(otra.coeficientes),
        )
        coeficientes = []
        for i in range(cogrado, indice_final):
            coef1 = self.coeficientes[i - self.cogrado] if self.cogrado <= i < self.cogrado + len(self.coeficientes) else 0
            coef2 = otra.coeficientes[i - otra.cogrado] if otra.cogrado <= i < otra.cogrado + len(otra.coeficientes) else 0
            coeficientes.append(coef1 + coef2)
        return SerieConPrincipio(coeficientes, cogrado)
    
    def __mul__(self, escalar_o_otra):
        """
        Producto por un escalar o producto convolución entre series.
    
        Parámetros:
        - escalar_o_otra (float | SerieConPrincipio): Escalar o otra serie.
    
        Retorna:
        - SerieConPrincipio: La serie resultante del producto.
        """
        if isinstance(escalar_o_otra, (int, float, Fraction)):
            coeficientes = [coef * escalar_o_otra for coef in self.coeficientes]
            return SerieConPrincipio(coeficientes, self.cogrado)
        elif isinstance(escalar_o_otra, SerieConPrincipio):
            cogrado = self.cogrado + escalar_o_otra.cogrado
            coeficientes = [0.0 if self.usa_float else Fraction(0)] * (len(self.coeficientes) + len(escalar_o_otra.coeficientes) - 1)
            for i, coef1 in enumerate(self.coeficientes):
                for j, coef2 in enumerate(escalar_o_otra.coeficientes):
                    coeficientes[i + j] += coef1 * coef2
            return SerieConPrincipio(coeficientes, cogrado)
        else:
            raise TypeError("El operador * solo admite un escalar o otra SerieConPrincipio.")
    
    def __rmul__(self, escalar):
        """
        Producto por un escalar desde la izquierda.
    
        Parámetros:
        - escalar (float): Escalar para multiplicar la serie.
    
        Retorna:
        - SerieConPrincipio: La serie resultante del producto.
        """
        return self.__mul__(escalar)
    def inversa(self, num_terminos=5):
        """
        Calcula la inversa de la serie.
    
        Parámetros:
        - num_terminos (int): Número de términos de la inversa a calcular.
    
        Retorna:
        - SerieConPrincipio: La serie inversa.
        """
        if not self.coeficientes or self.coeficientes[0] == 0:
            raise ValueError("La primera componente no nula debe ser distinta de cero.")
    
        b = [0.0 if self.usa_float else 0] * num_terminos
        a0 = self.coeficientes[0]
    
        if self.usa_float:
            b[0] = 1.0 / a0
        else:
            try:
                if isinstance(a0, (int, Fraction)):
                    b[0] = Fraction(1, a0)
                else:
                    b[0] = 1 / a0
            except Exception as e:
                raise ValueError(f"No se pudo calcular el inverso del término inicial {a0}: {e}")
    
        for j in range(1, num_terminos):
            suma = sum(
                b[r] * self.coeficientes[j - r]
                for r in range(j)
                if j - r < len(self.coeficientes)
            )
            b[j] = -b[0] * suma
       
        return SerieConPrincipio(b, -self.cogrado, final=True)
    
    def plot_serie(self, indices_previos_al_cogrado=3, title="Serie con principio (con cogrado)"):
        """
        Dibuja la serie como un gráfico de barras.
        """
        indices_nulos = range(self.cogrado - indices_previos_al_cogrado, self.cogrado)
        valores_nulos = [0] * indices_previos_al_cogrado  # Coeficientes nulos
    
        indices = range(self.cogrado, self.cogrado + len(self.coeficientes))
        valores = self.coeficientes
    
        # Combina valores nulos y de la serie
        all_indices = list(indices_nulos) + list(indices)
        all_valores = valores_nulos + valores
    
        fig = plt.figure(figsize=(9, 3))
        plt.stem(all_indices, all_valores)
        plt.axhline(0, color='black', linewidth=0.8)
    
        # Etiquetas para los índices
        etiquetas = [str(int(idx)) for idx in all_indices]
        plt.xticks(all_indices, etiquetas)  # Establecer etiquetas en el eje X
    
        # Personalizar las etiquetas
        for i, idx in enumerate(all_indices):
            xtick_label = plt.gca().get_xticklabels()[i]  # Obtener la etiqueta actual
            if idx == self.cogrado:
                xtick_label.set_color('red')  # Cambiar color a rojo
                xtick_label.set_weight('bold')  # Cambiar a negrita
            elif idx < self.cogrado:
                xtick_label.set_color('#D3D3D3')  # Cambiar color a gris claro
    
        plt.xlabel("Índice")
        plt.ylabel("Valor")
        plt.title(title)
        plt.tight_layout()
        plt.close(fig)
        return fig
        
    def subs(self, reglasDeSustitucion=[]):
        """ Sustitución de variables simbólicas """
    
        def CreaLista(t):
            """Devuelve t si t es una lista; si no devuelve la lista [t]"""
            return t if isinstance(t, list) else [t]
    
        def sustitucion(elemento, regla_de_sustitucion):
            return sympy.S(elemento).subs(CreaLista(regla_de_sustitucion))
            
        coeficientes = [sustitucion(elemento, reglasDeSustitucion) for elemento in self.coeficientes]
        return SerieConPrincipio(coeficientes, self.cogrado)    
    
class SerieConFinal:
    """
    Clase para representar series con un índice final arbitrario (grado).

    Atributos:
    - coeficientes (list): Lista de coeficientes sin ceros finales.
    - grado (int): Índice del último coeficiente no nulo.
    - usa_float (bool): Indica si la serie opera con números flotantes.
    - principio (bool): Indica si el principio de la serie se representa con una elipsis.

    Métodos:
    - __repr__: Representa la serie como un par (grado, lista de coeficientes).
    - _repr_latex_: Representa la serie en formato LaTeX.
    - __add__: Suma de series con final.
    - __mul__: Producto por escalar o convolución.
    - __rmul__: Producto por escalar desde la izquierda.
    - inversa: Calcula la inversa generando coeficientes hacia índices negativos.
    """

    def __init__(self, coeficientes, grado=0, principio=False, RepEscalera=False):
        # Eliminar ceros finales
        ultimos_no_nulos = next((i for i in reversed(range(len(coeficientes))) if coeficientes[i] != 0), None)
        if ultimos_no_nulos is None:  # Serie nula
            self.coeficientes = []
            self.grado = 0
        else:
            self.coeficientes = coeficientes[:ultimos_no_nulos+1]
            self.grado = grado - (len(coeficientes) - ultimos_no_nulos - 1)

        self.principio = principio
        self.usa_float = any(isinstance(c, float) for c in self.coeficientes if not hasattr(c, 'is_number'))
        self.coeficientes = [
            float(c) if self.usa_float and isinstance(c, (int, float)) else Fraction(c) if isinstance(c, (int, float)) else c
            for c in self.coeficientes
        ]
        self.RepEscalera=RepEscalera

    def __repr__(self):
        return f"({self.grado}, {self.coeficientes})"
   
    def latex(self):
        if not self.coeficientes:
            return "$0$"
        terminos = []
        for i, coef in enumerate(self.coeficientes):
            indice = self.grado - (len(self.coeficientes) - 1 - i)
            if coef != 0:
                coef_sympy = Rational(coef.numerator, coef.denominator) if isinstance(coef, Fraction) else coef
                coef_latex = r"\left(" + latex(coef_sympy) + r"\right)" if  isinstance(coef_sympy, Add) else latex(coef_sympy)
                if coef_sympy == -1 and indice != 0:
                    coef_latex = "-"
                if coef_sympy == 1 and indice != 0:
                    coef_latex = ""
                if indice == 0:
                    terminos.append(f"{coef_latex}")
                elif indice == 1:
                    terminos.append(f"{coef_latex}z")
                else:
                    terminos.append(f"{coef_latex}z^{{{indice}}}")
                    
        cadena = r'\cdots + ' if self.principio else ''
        representacion = cadena + " + ".join(terminos).replace("+ -", "- ")
        return representacion.replace("+ -", "- ")
        #return '$' + representacion.replace("+ -", "- ") + '$'
        #return '$'+ " + ".join(terminos).replace("+ -", "- ") + cadena +'$'


    def _repr_latex_(self):
        return '$'+self.latex()+'$'

    def _repr_html_(self):
        """ Construye la representación para el entorno jupyter notebook """
        """ Construye la representación para el entorno jupyter notebook """
        if self.RepEscalera:
            return html(self.escalera())
        else:
            return html(self.latex())
    
#    def _repr_png_(self):
#        """ Representación png para el entorno jupyter en Emacs """
#        try:
#            expr = '$'+self.latex()+'$'
#            workdir = tempfile.mkdtemp()
#            with open(join(workdir, 'borrame.png'), 'wb') as outputfile:
#                sympy.preview(expr, viewer='BytesIO', outputbuffer=outputfile)
#            return open(join(workdir, 'borrame.png'),'rb').read()
#        except:
#            return '$'+self.latex()+'$'

    def coeficientesEscalera(self):
        coeficientes = list(self.coeficientes)  # Creamos una copia de la lista de coeficientes
        
        if self.grado + 1 > len(self.coeficientes):
            coeficientes = [0]*(self.grado + 1 - len(self.coeficientes)) + coeficientes
            origen = 0            
        elif len(self.coeficientes) >  self.grado:
            coeficientes = coeficientes + [0]*(-self.grado)
            origen = len(coeficientes) - max(self.grado,0) -1 

        return coeficientes, origen
    
            
    def escalera(self):
        coeficientes, origen = self.coeficientesEscalera()
        celdas = [r'\cdots'] + [r'{\color{blue}{\{' + str(c) + '\}}}' if i == origen else f"{str(c)}" for i,c in enumerate(coeficientes)] + [r'0\cdots']        
        size = len(celdas)-1
        contenido = ' & '.join(celdas)
        latex_representation = "\\begin{array}{" + "c|"*size + "c}\n\\hline\n"
        latex_representation += f"{contenido} \\\\\hline"
        latex_representation += r"\end{array}"
        return '$' + latex_representation + '$'


    def __add__(self, otra):
        grado = max(self.grado, otra.grado)
        inicio = min(self.grado - len(self.coeficientes) + 1, otra.grado - len(otra.coeficientes) + 1)
        coeficientes = []
        for i in range(inicio, grado + 1):
            coef1 = self.coeficientes[i - (self.grado - len(self.coeficientes) + 1)] if self.grado - len(self.coeficientes) + 1 <= i <= self.grado else 0
            coef2 = otra.coeficientes[i - (otra.grado - len(otra.coeficientes) + 1)] if otra.grado - len(otra.coeficientes) + 1 <= i <= otra.grado else 0
            coeficientes.append(coef1 + coef2)
        return SerieConFinal(coeficientes, grado)

    def __mul__(self, escalar_o_otra):
        if isinstance(escalar_o_otra, (int, float, Fraction)):
            coeficientes = [coef * escalar_o_otra for coef in self.coeficientes]
            return SerieConFinal(coeficientes, self.grado)
        elif isinstance(escalar_o_otra, SerieConFinal):
            grado = self.grado + escalar_o_otra.grado
            coeficientes = [0.0 if self.usa_float else Fraction(0)] * (len(self.coeficientes) + len(escalar_o_otra.coeficientes) - 1)
            for i, coef1 in enumerate(self.coeficientes):
                for j, coef2 in enumerate(escalar_o_otra.coeficientes):
                    coeficientes[i + j] += coef1 * coef2
            return SerieConFinal(coeficientes, grado)
        else:
            raise TypeError("El operador * solo admite un escalar o otra SerieConFinal.")

    def __rmul__(self, escalar):
        return self.__mul__(escalar)

    
    def inversa(self, num_terminos=5):
        """
        Calcula la inversa de la serie como otra SerieConFinal.
        """
        if not self.coeficientes or self.coeficientes[-1] == 0:
            raise ValueError("La última componente no nula debe ser distinta de cero.")

        a = self.coeficientes[::-1]  # Trabajamos como si fuera una SerieConPrincipio
        b = [0.0 if self.usa_float else 0] * num_terminos
        a0 = a[0]

        if self.usa_float:
            b[0] = 1.0 / a0
        else:
            try:
                if isinstance(a0, (int, Fraction)):
                    b[0] = Fraction(1, a0)
                else:
                    b[0] = 1 / a0
            except Exception as e:
                raise ValueError(f"No se pudo calcular el inverso del coeficiente constante {a0}: {e}")

        for j in range(1, num_terminos):
            suma = sum(b[k] * a[j - k] for k in range(j) if j - k < len(a))
            b[j] = -b[0] * suma

        # El resultado es otra SerieConFinal con coeficientes en orden inverso
        return SerieConFinal(b[::-1], -self.grado, principio=True)

    
    def plot_serie(self, indices_posteriores_al_grado=3, title="Serie con final (con grado)"):
        """
        Dibuja la serie como un gráfico de barras, extendiéndola con ceros.
        """
        # Calcular los índices para el eje X
        longitud_lista = len(self.coeficientes)
        indices = []
        
        for i in range(longitud_lista):
            indice_mapeado = self.grado - (longitud_lista - 1 - i)
            indices.append(indice_mapeado)
        
        # Añadir índices nulos después de self.grado
        for i in range(1, indices_posteriores_al_grado + 1):
            indices.append(self.grado + i)
        
        # Crear la lista de valores
        valores = self.coeficientes + [0] * indices_posteriores_al_grado  # Añadir ceros
        
        fig = plt.figure(figsize=(9, 3))
        plt.stem(indices, valores)
        plt.axhline(0, color='black', linewidth=0.8)
        plt.xticks(indices, [str(int(idx)) for idx in indices])  # Etiquetas enteras
        plt.xlabel("Índice")
        plt.ylabel("Valor")
        plt.title(title)
        # Crear etiquetas
        etiquetas = [str(int(idx)) for idx in indices]
        plt.xticks(indices, etiquetas)  # Establecer etiquetas posiciones iniciales
        # Personalizar la etiqueta del índice self.grado
        for i, etiqueta in enumerate(etiquetas):
            if indices[i] == self.grado:
                etiq_actual, = plt.gca().get_xticklabels()[i],   # Obtener la etiqueta actual
                etiq_actual.set_color('red')  # Cambiar a azul
                etiq_actual.set_weight('bold')  # Cambiar a negrita
            elif indices[i] > self.grado:
                etiq_actual, = plt.gca().get_xticklabels()[i],   # Obtener la etiqueta actual
                etiq_actual.set_color('#D3D3D3')  # Cambiar a gris claro
                etiq_actual.set_fontstyle('italic')  # Cambiar a cursiva
                
        plt.tight_layout()
        plt.close(fig)
        return fig
        
    def subs(self, reglasDeSustitucion=[]):
        """ Sustitución de variables simbólicas """

        def CreaLista(t):
            """Devuelve t si t es una lista; si no devuelve la lista [t]"""
            return t if isinstance(t, list) else [t]

        def sustitucion(elemento, regla_de_sustitucion):
            return sympy.S(elemento).subs(CreaLista(regla_de_sustitucion))
            
        coeficientes = [sustitucion(elemento, reglasDeSustitucion) for elemento in self.coeficientes]
        return SerieConFinal(coeficientes, self.grado)    
    
def EscalerasConvolucion(secuencia1,secuencia2):
    serie1, origen1 = secuencia1.coeficientesEscalera()
    serie2, origen2 = secuencia2.coeficientesEscalera()
    
    negativos1 = serie1[:origen1]
    negativos2 = serie2[:origen2]
    negativos = max(len(negativos1), len(negativos2))
    lista1neg = [0]*(negativos-len(negativos1)) + negativos1 # añadimos ceros si no hay suficientes coeficientes a la izquierda del origen
    lista2neg = [0]*(negativos-len(negativos2)) + negativos2 # añadimos ceros si no hay suficientes coeficientes a la izquierda del origen
    cabeceraNeg = [r'z^{-' + str(i+1) +r'}' for i in reversed(range(len(lista1neg))) ]
    
    lista1origen = [r'{\color{blue}{\{' + str(serie1[origen1]) + '\}}}'] # destacamos el coeficiente situado en el origen
    lista2origen = [r'{\color{blue}{\{' + str(serie2[origen2]) + '\}}}'] # destacamos el coeficiente situado en el origen
    cabeceraOrigen = [r'{\color{blue}{z^0}}']
    
    positivos1 = serie1[origen1+1:]
    positivos2 = serie2[origen2+1:]
    positivos = max(len(positivos1), len(positivos2))
    lista1pos = positivos1 + [0]*(positivos-len(positivos1)) #[0 if i>len(lista1)-1 else lista1[i] for i in range(positivos)]
    lista2pos = positivos2 + [0]*(positivos-len(positivos2)) #=[0 if i>len(lista2)-1 else lista2[i] for i in range(positivos)]
    cabeceraPos = [r'z^{' + str(i+1) +r'}' for i,_ in (enumerate(lista1pos)) ]
    
    coeficientes1 = lista1neg + lista1origen + lista1pos
    coeficientes2 = lista2neg + lista2origen + lista2pos
    cabecera = cabeceraNeg + cabeceraOrigen + cabeceraPos
    
    # Consigue el tamaño para la alineación de las columnas
    size = len(coeficientes1)-1
    
    # Comenzamos a escribir la matriz en LaTeX
    #latex_matrix = "\\begin{array}{|" + "c|"*size + "}\n" #\n\\hline\n"
    latex_matrix = "\\begin{array}{" + "c|"*size + "c}\n" #\n\\hline\n"
    latex_matrix += " & ".join(str(x) for x in cabecera) + " \\\\\\hline\\hline\n"
    latex_matrix += " & ".join(str(x) for x in coeficientes1) + " \\\\\\hline\\hline\n"
    latex_matrix += " & ".join(str(x) for x in coeficientes2) + " \\\\\\hline\n"
    latex_matrix += "\\end{array}"

    # Mostramos la matriz en el notebook de Jupyter
    #display(Math(r'\['+latex_matrix+r'\]'))
    return Latex(r'\['+latex_matrix+r'\]')
