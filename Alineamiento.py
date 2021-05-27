import numpy as np
import math

mapeo = { 'A':0, 'C':1, 'G':2, 'T':3}

allineaciones = []

def get_caminos (cadena1, cadena2, pos_x, pos_y, matrix_similitud, matrix_matching, costo, grupos):
  if pos_x + pos_y == 0:
    allineaciones.append(grupos)
    return
  else:
    salida1 = matrix_matching[pos_y-1, pos_x-1] + matrix_similitud[mapeo[cadena1[pos_x-1]], mapeo[cadena2[pos_y-1]]] if pos_x-1 >=0 and pos_y-1 >=0 else -100000
    salida2 = matrix_matching[pos_y, pos_x-1] + costo if pos_x-1 >=0 else -100000
    salida3 = matrix_matching[pos_y-1, pos_x ] + costo if pos_y-1 >= 0 else -100000
    # Verificar posX y posY ya que cadenas depende de ello
    

    if salida1 == matrix_matching[pos_y, pos_x]:
      helper = grupos[:]
      helper[0] += cadena1[pos_x-1]
      helper[1] += cadena2[pos_y-1]
      get_caminos(cadena1, cadena2, pos_x-1, pos_y-1, matrix_similitud, matrix_matching, costo, helper)

    
    if salida2 == matrix_matching[pos_y, pos_x]:
      helper1 = grupos[:]
      helper1[0] += cadena1[pos_x-1]
      helper1[1] += "-"
      get_caminos(cadena1, cadena2, pos_x-1, pos_y, matrix_similitud, matrix_matching, costo, helper1)

    if salida3 == matrix_matching[pos_y, pos_x]:
      helper2 = grupos[:]
      helper2[0] += "-"
      helper2[1] += cadena2[pos_y-1]
      get_caminos(cadena1, cadena2, pos_x, pos_y-1, matrix_similitud, matrix_matching, costo, helper2)

def solution_global(cadena1, cadena2, matrix_similitud, costo):
  tam_cadena1 = len(cadena1)
  tam_cadena2 = len(cadena2)
  
  matriz_matching = np.zeros((tam_cadena2 + 1, tam_cadena1 + 1))
  for i in range(1, tam_cadena1+1):
    matriz_matching[0, i] = matriz_matching[0, i-1]  + costo

  for i in range(1, tam_cadena2+1):
    matriz_matching[i, 0] = matriz_matching[i-1, 0]  + costo

 
  for i in range(1, tam_cadena2 + 1):
    for j in range(1, tam_cadena1 + 1):
      matriz_matching[i,j] = max(matriz_matching[i-1, j-1] + matrix_similitud[mapeo[cadena1[j-1]], mapeo[cadena2[i-1]]], matriz_matching[i-1, j] + costo, matriz_matching[i, j-1] + costo)
  
  return matriz_matching

secuencia1 = input("Secuencia 1 ")
secuencia2 = input("Secuencia 2 ")
costo = input("Costo ")
costo = int(costo)

matrix_similitud = np.matrix(([[2, -7, -5, -7], [-7, 2, -7, -5], [-5, -7, 2, -7], [-7, -5, -7, 2]]))
print("Matriz de similitud")
print(matrix_similitud)

print("Matriz DP")
answer = solution_global(secuencia1, secuencia2, matrix_similitud, costo)
print(answer)

allineaciones = []

print("Alineaciones")
get_caminos (secuencia1, secuencia2, len(secuencia1), len(secuencia2), matrix_similitud, answer, costo, ["", ""])
for val in allineaciones:
  print(str(val[0][::-1]) + " " + str(val[1][::-1]))
print("Score")
print(answer[3][3])
