import OrbitaPotencialCorposMassivos
import OrbitaPotencialRaioDeLuz

tipo_orbita = input("Escolha 'M' para órbitas de corpos massivos e 'L' para órbitas de raios de luz: ")

if tipo_orbita == "M":

    valor = input("Insira o valor do momento angular adimensional: ")
    orbitaCorposMassivos = OrbitaPotencialCorposMassivos.OrbitaPotencialCorposMassivos(valor)

    if (orbitaCorposMassivos.calcularCorposMassivos()):
        print("Calculando corpos massivos")
        orbitaCorposMassivos.calcularPlotPotencialVisaoUm()
        energia = input("Insira o valor do parâmetro de energia: ")
        numeroOrbitas = input("Para uma órbita ligada (1 ≤ E < 20), escolha também o número de órbitas que deseja traçar: ")
        orbitaCorposMassivos.calcularPlotPotencialVisaoDois(energia, numeroOrbitas)
    else:
        print("Para esse valor do momento angular, a energia potencial efetiva não possui um mínimo ou máximo local.")
        orbitaCorposMassivos.calcularPlotPotencialVisaoTres()
        energia = input("Insira o valor do parâmetro de energia: ")
        numeroOrbitas = input("Para uma órbita ligada (1 ≤ E < 20), escolha também o número de órbitas que deseja traçar: ")
        orbitaCorposMassivos.calcularPlotPotencialVisaoQuatro(energia, numeroOrbitas)

else:
    impacto = input("Escolha o valor do parâmetro de impacto 'd', em unidades de $r_g$: ")
    orbitaPotencialRaioDeLuz = OrbitaPotencialRaioDeLuz.OrbitaPotencialRaioDeLuz(impacto)

    orbitaPotencialRaioDeLuz.calcular_orbita_potencial_raio_luz()
        