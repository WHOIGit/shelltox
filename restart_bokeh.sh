killall bokeh

nohup /opt/miniconda3/envs/shelltox/bin/bokeh serve bokeh_app/shelltox.py --address westerly.whoi.edu --port 5100 --allow-websocket-origin=westerly.whoi.edu:5100 2>&1 > bokeh_app/logs/bokeh_serve.log &
