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


Files:
- stats.xlsx: Este archivo contiene diferentes datos sobre el generador y las rutas. Tenemos el número de metabolitos que hay en cada step del generador, cuantos de estos ya han pasado y cuantos faltan por pasar. También tiene una lista de los matches con E.Coli y producibles. Sobre las rutas, tiene la cantidad de pairs en cada diámetro, los metabolitos del grafo y los detectables y producibles que aparecen en el grafo. Hay dos versiones de esta tabla, porque haciendo pruebas he encotrado otra manera de calcular el grafo que obtiene muchas mas pairs.
