PYUI=pyuic5
NAME_LAYOUT=QtLayout
NAME_PARAM=InputParameter

all : layout

.PHONY: all 

layout:
	$(PYUI) $(NAME_LAYOUT).ui -o $(NAME_LAYOUT).py

input_param:
	$(PYUI) $(NAME_PARAM).ui -o $(NAME_PARAM).py

clean:
	rm -f $(NAME).py
