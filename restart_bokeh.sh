killall bokeh

nohup bokeh serve bokeh_app/model_app.py --address shelltox.whoi.edu --port 5100 --allow-websocket-origin=shelltox.whoi.edu 2>&1 > bokeh_app/logs/bokeh_serve.log &

