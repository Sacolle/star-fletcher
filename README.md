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

Por padrão o programa usa um limite de 200000 tarefas pendentes (ao passar disso, a submissão
espera até restarem 180000), definido em `DEFAULT_LIMIT_{MAX,MIN}_SUBMITTED_TASKS` no `main.c`.
Esses valores são exportados como `STARPU_LIMIT_{MAX,MIN}_SUBMITTED_TASKS` antes do `starpu_init`,
mas só quando a variável ainda não existe no ambiente, então o que o usuário exportar tem prioridade.
Com ele o 200³ (8 CPUs) fica em ~3.6 GB; sem ele passa de 24 GB e pode travar a máquina. Esse limite
também define quanto overlap entre iterações é possível, então é um parâmetro de experimento.
`STARPU_LIMIT_MAX_SUBMITTED_TASKS=-1` desliga o limite. O StarPU só bloqueia quando o número de tarefas
pendentes passa dos dois valores, então exportar só um `MAX` menor que o `MIN` padrão faz o `MIN`
valer como limite; o ideal é definir os dois:

```bash
STARPU_LIMIT_MAX_SUBMITTED_TASKS=50000 STARPU_LIMIT_MIN_SUBMITTED_TASKS=45000 ./main ...
```

Outra possibilidade é utilizar `starpu_task_wait_for_n_submitted()` para que o loop de execução
bloqueie depois de submeter `n` tarefas, limitando em iterações em vez de tarefas.


### GPU (CUDA)

A partição dos cubos precisa do StarPU master (no 1.4.12 o resultado sai errado), então a versão
CUDA é compilada no shell `cuda-lattest`:

```bash
nix develop .#cuda-lattest
make CUDA_BACKEND=1 RELEASE_MODE=1                  # kernel na CPU e na GPU
make CUDA_BACKEND=1 RELEASE_MODE=1 NO_CPU_KERNEL=1  # tarefas RTM só na GPU
```

`ARCH` escolhe o `-arch` do nvcc (ex. `ARCH=sm_80`). Sem ele, o Makefile usa a compute capability
da GPU visível (`nvidia-smi`) e para com erro se não achar nenhuma, em vez de compilar para uma
arquitetura antiga como faz o `-arch=native` do nvcc. O kernel precisa de sm_70 ou mais novo. Para
o `nix build` dos pacotes CUDA (sem GPU no sandbox), passe `cuda_arch = "sm_XX"` no `star-fletcher.nix`.

Dentro do `nix develop` numa máquina que não é NixOS, o driver da GPU do sistema (`libcuda.so`) não
fica visível: o StarPU não encontra a GPU (as tarefas RTM falham com `No such device` quando o kernel
de CPU está desligado). Rode o programa pelo `nixglhost`, incluído no shell `cuda-lattest`; no NixOS,
use `LD_LIBRARY_PATH=/run/opengl-driver/lib`:

```bash
STARPU_NCPU=1 STARPU_NCUDA=1 nixglhost ./main TTI ...
```

As tarefas de perturbação e de escrita só existem na CPU, então rode com pelo menos um worker de CPU.

`EXACT_FP=1` compila sem fused multiply-add na CPU (`-ffp-contract=off`) e na GPU (`-fmad=false`),
para que as duas façam as mesmas operações e a saída seja bit a bit igual. Os MD5 de referência
foram gravados em x86-64, onde o gcc não funde operações; em aarch64 ele funde por padrão e a libm é
outra, então lá os MD5 são diferentes. Por isso as comparações são feitas na mesma máquina.

`scripts/validate-gpu.sh` compara, na mesma máquina, a CPU (`EXACT_FP=1`) com a GPU: MD5 idêntico no
build exato e diferença máxima frame a frame no build normal.

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
