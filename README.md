# SME Robotics Engine

**Autor**: Luis Santos aka DJOKER  
**Versão**: 0.1  

Um motor de simulação física e visualização 3D em tempo real usando **Raylib**, **ODE (Open Dynamics Engine)** e **Pybind11** para ligação com **Python**.

---

## 🎯 Objetivo

Este projeto permite criar simulações físicas interativas com rendering em 3D, oferecendo:

- Física realista com ODE
- Visualização com Raylib
- Prototipagem rápida com Python
- Hot reload de scripts Python em tempo real
- Criação de rigid bodies, veículos e juntas físicas via bindings

---



## 🔧 Requisitos

- C++17
- Python ≥ 3.6
- Raylib
- ODE (Open Dynamics Engine)
- Pybind11

```bash
sudo apt install libode-dev libpython3-dev python3-pybind11
```

---

## ⚙️ Compilação

```bash
mkdir build
cd build
cmake ..
make
```

---

## 🚀 Execução

```bash
./sme ../main.py
```

---

## 🧪 Exemplo Python: `main.py`

```python
from sme import Cube, SetGravity

WINDOW_WIDTH = 1280
WINDOW_HEIGHT = 720
WINDOW_TITLE = "SME Teste"
HOT_RELOAD = True

cube = None

def load():
    global cube
    SetGravity(0, -9.81, 0)
    cube = Cube(1.0, 1.0, 1.0, 1.0)

def unload():
    print("Unload OK")

def loop():
    pass
```

---

## 🧩 Bindings Disponíveis

### Objetos Físicos:

- `Cube`, `Sphere`, `Cylinder`, `Capsule`
- `Vehicle` com `.AddWheel()`, `.SetSuspension()`, `.SetSteering()`
- Juntas: `HingeJoint`, `BallJoint`, `FixedJoint`, `SliderJoint`, `Hinge2Joint`, `MotorJoint`

### Utilidades:

- `SetGravity(x, y, z)`
- `ClearWorld()`
- Inputs: `IsKeyPressed(KEY_W)`, `GetMousePosition()`
- Rendering: `DrawFPS(x, y)`

---

## 🌀 Hot Reload

Sempre que alterares o `main.py`, o motor deteta automaticamente e recarrega o script sem reiniciar a aplicação.

---

## 📸 Screenshot

> (Adiciona aqui uma imagem da simulação em execução)

---

## 📜 Licença

Uso livre para fins educacionais e prototipagem experimental.

---

## 📣 Desenvolvido por

Luis Santos  
42 Porto / SEA.ME  
2025
