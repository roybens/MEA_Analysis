from dash import Dash
from layout import create_layout
from app_callbacks import register_callbacks

app = Dash(__name__, suppress_callback_exceptions=True)
app.layout = create_layout()

register_callbacks(app)

def find_free_port():
    import socket
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.bind(('', 0))
    port = s.getsockname()[1]
    s.close()
    return port 

if __name__ == '__main__':
    free_port = find_free_port()
    app.run_server(debug=True, port=free_port)
