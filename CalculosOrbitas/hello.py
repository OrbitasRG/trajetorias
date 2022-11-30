from flask import Flask, request, jsonify
from flasgger import Swagger
from flask_restful import Api, Resource

import OrbitaCeleste
import OrbitaRaioDeLuz

app = Flask(__name__)
swagger = Swagger(app)

@app.route('/orbita/raio_de_luz/<int:posicaoInicial>/<int:velocidadeInicial>/')
def raio_de_luz(posicaoInicial, velocidadeInicial):
    """Calcula o raio de luz
    Faz o cálculo do raio de luz e retorna um gif com a operação realizada
    ---
    parameters:
      - name: posicaoInicial
        in: path
        type: int
        required: false
      - name: velocidadeInicial
        in: path
        type: int
        required: false
    responses:
      200:
        description: Retorna um string png com a animação
        schema:
        examples:
    """
    print("Entrou para calcular")
    orbitaCeleste = OrbitaCeleste.OrbitaCeleste()
    return orbitaCeleste.calcular_celeste(posicaoInicial, velocidadeInicial)

app.run(debug=True)