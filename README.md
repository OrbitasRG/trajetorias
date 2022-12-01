# Trajetorias

Este projeto visa mover o streamlit para um servidor com mais potência

# Instalação

Para criar as imagens, vá até as pastas GIFs e Potencial

Na pasta GIFs execute o seguinte comando

docker build -t gifs .

docker run -p 8502:8502 gifs

Na pasta Potencial execute o seguinte comando

docker build -t potencial .

docker run -p 8501:8501 potencial

Ele irá criar duas imagens, com elas conseguimos subir em um servidor VPS ou até local para fazer a execução do projeto