# Star-Fletcher

Star-Fletcher é uma implementação da propagação de onda com o tempo do algoritmo de Reverse Time Migration (RTM) [Fletcher](https://doi.org/10.1190/1.3269902).
Ele utiliza paralelismo de tarefas em primeira ordem com o framework [StarPU](https://starpu.gitlabpages.inria.fr/). 
O flake [nix-starpu](https://github.com/Sacolle/nix-starpu) é usado para compilar esse programa em sistemas nix.

Utilizando o [Madagascar](https://github.com/Sacolle/nix-madagascar) para a visualização do meio e o [Eztrace](https://github.com/Sacolle/eztrace-pallas-nix) para rastros 

## Execução

Este é um projeto nix, então se tiver o package manager instalado ou usar o OS, 
basta rodar `nix develop` para entrar no ambiente de desenvolvimento. Após isso,
`make` compila o código e `make run` roda com uma entrada de exemplo.


### Limite de submissão de tarefas

Devido ao uso de partições, cada cubo no volume agora constitui de 19 `data_handle_t`, 
um original e 18 das partições. Dependendo da segmentação do espaço, esses `data_handle_t`
podem acumular e consumir bastante memória. A medida que o programa executa, eles são 
descartados, então podemos colocar um limite no número de tarefas submetidas, assim, não
sobrecarregamos a memória.

Por enquanto o limite é passado pelo ambiente:

```bash
STARPU_LIMIT_MAX_SUBMITTED_TASKS=200000 STARPU_LIMIT_MIN_SUBMITTED_TASKS=180000 ./main ...
```

Com isso o 200³ (8 CPUs) fica em ~3.6 GB. Esse limite também define quanto overlap entre iterações
é possível, então é um parâmetro de experimento.

Há a possibilidade de (a) adicionar como configuração na chamada de `starpu_init` e de (b) utilizar
`starpu_task_wait_for_n_submitted()` para que o loop de execução bloqueie depois de submeter `n` tarefas.


## Uso no Emacs

Agora foi adicionado dois arquivos,
`compile_commands.json` and `.envrc`.
O primeiro é usado pelo `eglot` para fazer todo o paranauê do lsp funcionar.
O segundo serve para o emacs automaticamente entrar no dev enviroment quando eu estou acessando um arquivo naquele projeto.


Caso adicione um arquivo, use
```bash
make lsp
```

e rode `M-x eglot-reconnect` para ter certeza que o LSP atualizou
