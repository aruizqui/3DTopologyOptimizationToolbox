# transformData.py
# Transform Optistruct '.fem' file to coordinates and 
# node-Connectivity text files. Requires filename argument

import sys
import os


def transform_data(file_path):

    # Obtain root where the python script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Output files in same file as the script
    coords_file= os.path.join(script_dir,'coordinates.txt')
    connec_file= os.path.join(script_dir,'conectivity.txt')
    
    # Crea los archivos de salida para coordenadas y conectividad
    with open(file_path, 'r') as file:
        #Create coordinates and conectivity files.
        with open(coords_file,'w') as fcoord, open(connec_file, 'w') as fconn:

            #Temp list for HEX8 nodes
            chexa_nodes = [] 

            for line in file: 
            #Procesa elementos GRID
              if line.startswith("GRID"):
                x2 = str(float(line[16:32]))
                x3 = str(float(line[32:40]))
                x4 = str(float(line[40:49]))
                fcoord.write(f"{x2} {x3} {x4}\n")

            #Procesa elementos CTETRA
              if line.startswith("CTETRA"):
                x2 = str(int(line[25:33]))
                x3 = str(int(line[33:41]))
                x4 = str(int(line[41:49]))
                x5 = str(int(line[49:57]))
                fconn.write(f"{x2} {x3} {x4} {x5}\n")

            #Procesa elementos CHEXA  
              if line.startswith("CHEXA"):
                chexa_nodes = [
                  str(int(line[25:33])),
                  str(int(line[33:41])),
                  str(int(line[41:49])),
                  str(int(line[49:57])),
                  str(int(line[58:65])),
                  str(int(line[65:72])),
                ]
              elif line.startswith('+') and chexa_nodes:
                chexa_nodes += [
                  str(int(line[8:17])),
                  str(int(line[17:25]))
                  ]
              if len(chexa_nodes)== 8:
                  fconn.write(" ".join(chexa_nodes) + "\n")
                  chexa_nodes = []

if __name__ == "__main__":
    if len(sys.argv) !=2:
        print("Error: Se debe proporcionar la ruta al archivo .bdf")
        sys.exit(1)
    
    file_path = sys.argv[1]
    transform_data(file_path)
