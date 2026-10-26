"""Build and run an isolated v1706 case; never modify the tracked case."""
import argparse
from pathlib import Path
import re
import shutil
import subprocess

root=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser()
parser.add_argument('--startup',action='store_true')
parser.add_argument('--mesh-only',action='store_true')
args=parser.parse_args()
work=root/'.work'/('startup' if args.startup else 'baseline')
if work.exists(): raise SystemExit('Output already exists: '+str(work))
work.mkdir(parents=True)
for name in ['0','constant','system']: shutil.copytree(root/'case'/name,work/name)
shutil.rmtree(work/'constant/polyMesh')
def run(argv,name):
    with (work/name).open('w') as log:
        subprocess.run(argv,cwd=work,stdout=log,stderr=subprocess.STDOUT,check=True)
run(['gmsh',str(root/'geometry/su2Wedge.geo'),'-3','-format','msh2','-o','wedge.msh'],'log.gmsh')
run(['gmshToFoam','wedge.msh'],'log.gmshToFoam')
p=work/'constant/polyMesh/boundary';s=p.read_text()
for name,kind in [('inlet','patch'),('outlet','patch'),('upper','wall'),('lower','wall'),('frontAndBack','empty')]:
    s,count=re.subn(r'(\b'+name+r'\s*\{[^}]*?\btype\s+)\w+(;)',r'\g<1>'+kind+r'\2',s,flags=re.S)
    assert count==1,name
p.write_text(s)
if args.startup:
    p=work/'system/controlDict';s=p.read_text()
    s=re.sub(r'endTime\s+[^;]+;','endTime 0.01;',s)
    s=re.sub(r'writeInterval\s+[^;]+;','writeInterval 0.00005;',s)
    p.write_text(s)
(work/'wedge.foam').touch()
run(['checkMesh'],'log.checkMesh')
if 'Mesh OK.' not in (work/'log.checkMesh').read_text(): raise SystemExit('Mesh check did not pass')
if not args.mesh_only: run(['rhoCentralFoam'],'log.rhoCentralFoam')
print(work)
