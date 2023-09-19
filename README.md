# metabolic-map

Scripts:
- generator.py: Script para generar los nuevos metabolitos
- smile_normalizer.py: Script que normaliza los smiles en el formato común a todo el proyecto
- get_all_participants.py: Script que obtiene todos los participantes del grafo a partir de los resultados obtenidos del generador
- create_graph.py: Script que a partir de una lista de metabolitos crea el grafo que los une
- find_path.py: Script con las funciones para obtener todos los pathways entre dos grupos de metabolitos
- expand_path.py: Script con las funciones para obtener todos los pathways con sus diferentes rutas y bifurcaciones
- path_clenaer.py: Script para eliminar posibles paths duplicados dentro de un mismo pathways
- scope.py: Script con las funciones para obtener el scope de una ruta
