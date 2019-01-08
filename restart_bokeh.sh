killall bokeh

nohup bokeh serve model_app.py --address shelltox.whoi.edu --port 5100 --allow-websocket-origin=shelltox.whoi.edu:5100 --allow-websocket-origin=shelltox.whoi.edu 2>&1 > logs/bokeh_serve.log &

