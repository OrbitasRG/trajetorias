import OrbitaCeleste
import OrbitaRaioDeLuz

tipo_orbita = input("Escolha 'M' para órbita de corpos celestes ou 'L' órbitas de raios de luz: ")

if tipo_orbita == "M":
    x0 = input("Escolha o valor da posição inicial (em km): ")
    v0 = input("Escolha o valor da velocidade inicial (em unidades da velocidade da luz): ")

    orbitaCeleste = OrbitaCeleste.OrbitaCeleste()
    orbitaCeleste.calcular_celeste(x0, v0)

if tipo_orbita == "L":
    v0 = input("Escolha o valor do parâmetro de impacto $d$ (em km): ")

    orbitaRaioDeLuz = OrbitaRaioDeLuz.OrbitaRaioDeLuz()
    orbitaRaioDeLuz.calcular(v0)

