"""Encode the analysis PNG frames with the available GStreamer VP9 encoder."""
from pathlib import Path
import subprocess
root=Path(__file__).resolve().parents[1]
frames=sorted((root/'media/frames').glob('frame-*.png'))
if not frames: raise SystemExit('No rendered frames')
assert [p.name for p in frames]==[f'frame-{i:04d}.png' for i in range(len(frames))]
subprocess.run(['gst-launch-1.0','-e','multifilesrc','location='+str(root/'media/frames/frame-%04d.png'),'index=0','stop-index='+str(len(frames)-1),'caps=image/png,framerate=20/1','!','pngdec','!','videoconvert','!','vp9enc','deadline=1','cpu-used=6','target-bitrate=900000','!','webmmux','!','filesink','location='+str(root/'media/startup-plots.webm')],check=True)
