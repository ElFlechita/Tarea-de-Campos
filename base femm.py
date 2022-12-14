# Se Importan las funciones a usar de las librerias requeridas

from femm import  ei_addarc, ei_addblocklabel, ei_addboundprop, eo_makeplot, eo_gete, ei_addmaterial, ei_addnode, eo_savebitmap, ei_addsegment, ei_analyze, ei_clearselected, ei_loadsolution, ei_probdef, ei_saveas, ei_selectarcsegment, ei_selectlabel, ei_setarcsegmentprop, ei_setblockprop, eo_showcontourplot, eo_zoom, hideconsole, newdocument, openfemm, eo_hidepoints, eo_maximize, eo_hidemesh, eo_hidegrid, eo_showdensityplot
from numpy import linspace, concatenate
from numpy import abs as ABS
from math import pow

# Se define la clase Permitividad, la cual sirve para crear las distintas permitividades en función del radio.
#
# Forma de uso:
#
# obj = Permitividad()
#
# En caso de necesitar un perfil de permitividad lineal, usar:
#
# f_lineal = obj.lineal(r_min = 0.01, r_max = 0.015, epsilon_min = 0.2, epsilon_max = 5)
#
# En caso de necesitar un perfil de permitividad cuadrática, usar:
#
# f_cuadratica = obj.cuadratica(r_min = 0.01, r_max = 0.015, epsilon_min = 0.2, epsilon_max = 5)
#
# En caso de necesitar un perfil de permitividad asintótico (fracción), usar:
#
# f_fraccion = obj.fraccion(r_min = 0.01, r_max = 0.015, epsilon_min = 0.2, epsilon_max = 5)
#
# En caso de necesitar un perfil de permitividad constante, usar:
#
# f_constante = obj.constante(const = 1, r_min = 0.01, r_max = 0.015, epsilon_min = 0.2, epsilon_max = 5)
#
# Donde r_min, r_max, epsilon_min, epsilon_max y const son argumentos opcionales y estan definidos por defecto en las variables estándar del problema.

class Permitividad():
    def __init__(self):
        pass
    
    def cuadratica(self, r_min = 0.01, r_max = 0.015, epsilon_min = 0.2, epsilon_max = 5):
        return lambda r: (pow(r_max, 2) + pow(r, 2) * (-1 + epsilon_max) - pow(r_min, 2) * epsilon_max) / (pow(r_max, 2) - pow(r_min, 2))
    
    def lineal(self, r_min = 0.01, r_max = 0.015, epsilon_min = 0.2, epsilon_max = 5):
        return lambda r: (r_max + r * ( -1 + epsilon_max ) - r_min * epsilon_max) / (r_max - r_min)
    
    def fraccion(self, r_min = 0.01, r_max = 0.015, epsilon_min = 0.2, epsilon_max = 5):
        return lambda r: (-r * r_min + r_max * r_min + r * r_max * epsilon_min - r_max * r_min * epsilon_min)/(r * ( r_max - r_min ))
    
    def constante(self, const = 1, r_min = 0.01, r_max = 0.015, epsilon_min = 0.2, epsilon_max = 5):
        return lambda r: const
    
# Se define la función femm, la cual ejecuta la aplicación, y define los parametros del problema como argumentos modificables
#
# Forma de uso:
#
# femm( f1, file = "file", rmin = 0.01, rmax = 0.015, rdip = 0.0035, emax = 5, emin = 0.2, V = 0, n = 20, qs = 0.00001,
# unidad = "centimeters", simetria = 'axi', save_bitmap = False, get_e = False, generate_density_plot = (False, 1, 0, 0, 5, -5),
# generate_contour_plot = (False, 20, -25, 25) )

def femm(f1, file = "file", rmin = 0.01, rmax = 0.015, rdip = 0.0035, emax = 5, emin = 0.2, V = 0, n = 20, qs = 0.00001, unidad = "centimeters", simetria = 'axi', save_bitmap = False, get_e = False, generate_density_plot = (False, 1, 0, 0, 5, -5), generate_contour_plot = (False, 20, -25, 25)):
    openfemm()
    newdocument(1)
    ei_probdef(unidad, simetria, 1e-8,0,30)
    rmin, rmax, rdip, emax, emin, V, n, qs = rmin, rmax, rdip, emax, emin, V, n, qs
    puntos = concatenate((linspace(-rmax, -rmin, n+1),linspace(rmin, rmax, n+1)))
    paso = puntos[1] - puntos[0]
    for i in puntos:
        ei_addnode(0, i)
        if ABS(i) > rmin:
            if i > 0:
                ei_addmaterial(str(round(f1(i),2)), f1(i), f1(i), 0)
                ei_addblocklabel(0,i-paso/2)
                ei_selectlabel(0,i-paso/2)
                ei_setblockprop(str(round(f1(i),2)), 1, 0, 0)
                ei_clearselected()
    ei_addboundprop('dip_p', 0, qs, 0, 0, 2) 
    ei_addboundprop('dip_n', 0, -qs, 0, 0, 2) 
    ei_addboundprop('0v', V, 0, 0, 0, 0)
    for i in range(n+1):
        ei_addarc(0,puntos[i],0,puntos[-(1+i)],180,1)
        if i == 0:
            ei_selectarcsegment(0,puntos[i])
            ei_setarcsegmentprop(1, '0v', 0, 0, '<None>')
            ei_clearselected()
    ei_addmaterial('Vacio', 1, 1, 0)
    ei_addblocklabel(rmin/2,0)
    ei_selectlabel(rmin/2,0)
    ei_setblockprop('Vacio', 1, 0, 0)
    ei_clearselected()
    ei_addnode(0, rmin/2 + rdip)
    ei_addnode(0, rmin/2 - rdip)
    ei_addnode(0, -rmin/2 + rdip)
    ei_addnode(0, -rmin/2 - rdip)
    ei_addsegment(0, rmin/2 + rdip, 0, rmax)
    ei_addsegment(0, -rmin/2 - rdip, 0, -rmax)
    ei_addsegment(0, -rmin/2 + rdip, 0, rmin/2 - rdip)
    ei_addarc(0,rmin/2 - rdip,0,rmin/2 + rdip,180,1)
    ei_addarc(0,-rmin/2 - rdip,0,-rmin/2 + rdip,180,1) 
    ei_selectarcsegment(0,rmin/2)
    ei_setarcsegmentprop(1, 'dip_p', 0, 0, '<None>')
    ei_clearselected()
    ei_selectarcsegment(0,-rmin/2)
    ei_setarcsegmentprop(1, 'dip_n', 0, 0, '<None>')
    ei_clearselected()
    ei_saveas(f"{file}.FEE")
    ei_analyze(1)
    ei_loadsolution()
    eo_maximize()
    if generate_density_plot[0]:
        eo_showdensityplot(generate_density_plot[1], generate_density_plot[2], generate_density_plot[3], generate_density_plot[4], generate_density_plot[5])
    if generate_contour_plot[0]:
        eo_showcontourplot(generate_contour_plot[1], generate_contour_plot[2], generate_contour_plot[3])
    if save_bitmap:
        eo_savebitmap(f"{file[:-4]}.BMP")
    if get_e:
        ex0, ey0 = eo_gete(0.0001, 0.0098)
        ex1, ey1 = eo_gete(0.0005, 0.0101)
        return (ex0, ex1)    

# Se define el objeto Permitividad.

f = Permitividad()

# Se definen los perfiles estándar en el diccionario Permitividades.

Permitividades = {"Constante": f.constante(1),
                  "Fraccion": f.fraccion(),
                  "Lineal": f.lineal(),
                  "Cuadratico": f.cuadratica()}

# Se ejecutan los casos estándar.

femm(Permitividades["Lineal"], file = f"Caso lineal estandar")
femm(Permitividades["Cuadratico"], file = f"Caso cuadratico estandar")
femm(Permitividades["Fraccion"], file = f"Caso asintotico estandar")
femm(Permitividades["Constante"], file = f"Caso sin dieléctrico estandar")

# En el caso 3 se ejecutan iteraciones para el caso lineal y para el caso cuadrático, hasta que el parámetro val alcance un valor de 0.1

# Se definen listas vacias para usar en análisis.

test_lineal = []
test_cuadratica = []
cases_lineal = []
cases_cuadratica = []

for i in range(5, 95, 5):
    test_lineal.append(f.lineal(r_min = 0.01, r_max = 0.015, epsilon_min = 0.2, epsilon_max = i))
    cases_lineal.append(i)
    
for i in range(100, 1000, 100):
    test_lineal.append(f.lineal(r_min = 0.01, r_max = 0.015, epsilon_min = 0.2, epsilon_max = i))
    cases_lineal.append(i)
    
for i in range(1000, 4000, 1000):
    test_lineal.append(f.lineal(r_min = 0.01, r_max = 0.015, epsilon_min = 0.2, epsilon_max = i))
    cases_lineal.append(i)
    
for i in range(5000, 80000, 5000):
    test_lineal.append(f.lineal(r_min = 0.01, r_max = 0.015, epsilon_min = 0.2, epsilon_max = i))
    cases_lineal.append(i)
    
for i in range(80000, 82000, 200):
    test_lineal.append(f.lineal(r_min = 0.01, r_max = 0.015, epsilon_min = 0.2, epsilon_max = i))
    cases_lineal.append(i)
    
for i in range(82010, 82080, 10):
    test_lineal.append(f.lineal(r_min = 0.01, r_max = 0.015, epsilon_min = 0.2, epsilon_max = i))
    cases_lineal.append(i)
    
for i in range(82080, 82085, 1):
    test_lineal.append(f.lineal(r_min = 0.01, r_max = 0.015, epsilon_min = 0.2, epsilon_max = i))
    cases_lineal.append(i)
    
for i in range(5, 95, 5):
    test_cuadratica.append(f.cuadratica(r_min = 0.01, r_max = 0.015, epsilon_min = 0.2, epsilon_max = i))
    cases_cuadratica.append(i)
    
for i in range(100, 1000, 100):
    test_cuadratica.append(f.cuadratica(r_min = 0.01, r_max = 0.015, epsilon_min = 0.2, epsilon_max = i))
    cases_cuadratica.append(i)
    
for i in range(1000, 4000, 1000):
    test_cuadratica.append(f.cuadratica(r_min = 0.01, r_max = 0.015, epsilon_min = 0.2, epsilon_max = i))
    cases_cuadratica.append(i)
    
for i in range(5000, 100000, 5000):
    test_cuadratica.append(f.cuadratica(r_min = 0.01, r_max = 0.015, epsilon_min = 0.2, epsilon_max = i))
    cases_cuadratica.append(i)
    
for i in range(100000, 100600, 100):
    test_cuadratica.append(f.cuadratica(r_min = 0.01, r_max = 0.015, epsilon_min = 0.2, epsilon_max = i))
    cases_cuadratica.append(i)
    
for i in range(100600, 100640, 10):
    test_cuadratica.append(f.cuadratica(r_min = 0.01, r_max = 0.015, epsilon_min = 0.2, epsilon_max = i))
    cases_cuadratica.append(i)
    
for i in range(100640, 100650, 1):
    test_cuadratica.append(f.cuadratica(r_min = 0.01, r_max = 0.015, epsilon_min = 0.2, epsilon_max = i))
    cases_cuadratica.append(i)

def test_lineal():
    count = 0
    for test in test_lineal:
        ex0, ex1 = femm(test, file = f"TestLineal{count}")
        epsilon_i = cases_lineal[count]
        val = 100 * ex1 / ex0
        print(f"{epsilon_i},{val}")
        count += 1

def test_cuadratica():
    count = 0
    for test in test_cuadratica:
        ex0, ex1 = femm(test, file = f"TestCuadratico{count}")
        epsilon_i = cases_lineal[count]
        val = 100 * ex1 / ex0
        print(f"{epsilon_i},{val}")
        count += 1

test_lineal()
test_cuadratica()   


'''
for i in range(10, 1000, 10):
    qs = i / 1000000
    femm(Permitividades["Constante"], file = f"TestLineal{1}", qs = qs)
'''
    
#femm(f.lineal(r_min = 0.01, r_max = 0.015, epsilon_min = 0.2, epsilon_max = 82000), file = f"TestCuadratico{00}.FEE")
#femm(f.cuadratica(r_min = 0.01, r_max = 0.015, epsilon_min = 0.2, epsilon_max = 100647.3250), file = f"TestCuadratico{00}.FEE")


