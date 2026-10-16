{
    description = "A very basic flake";

    inputs = {
        nixpkgs.url = "github:nixos/nixpkgs?ref=nixos-26.05";
        StarPU = {
          url = "github:Sacolle/nix-starpu";
          inputs.nixpkgs.follows = "nixpkgs";
        };
        nix-gl-host = {
            url = "github:numtide/nix-gl-host";
            inputs.nixpkgs.follows = "nixpkgs";
        };
        eztrace = {
            url = "github:Sacolle/eztrace-pallas-nix";
            inputs.nixpkgs.follows = "nixpkgs";
        };

        #madagascar.url = "github:Sacolle/nix-madagascar";
        flake-utils.url = "github:numtide/flake-utils";
    };

    outputs = { self, nixpkgs, StarPU, eztrace, /* madagascar,*/ nix-gl-host, flake-utils }: 
    flake-utils.lib.eachSystem [ "x86_64-linux" "aarch64-linux" ] (system:
    let 
        pkgsconfigs = { 
            inherit system; 
            config.allowUnfree = true;
        };
        pkgs = import nixpkgs pkgsconfigs;

        cudaPacks = pkgs.cudaPackages_13_0;
        # the host compiler supported by this CUDA version (nvcc rejects newer gcc's)
        cudaStdenv = cudaPacks.backendStdenv;

        # StarPU 1.4.12 uses cudaDeviceProp::computeMode and two cuFFT error codes, all removed in CUDA 13.
        # starpu.nix overrides patchPhase, so `patches` would be ignored: apply it there instead
        starpu-base = StarPU.packages.${system}.default.overrideAttrs (old: {
            patchPhase = "patch -p1 < ${./patches/starpu-cuda13.patch}\n" + old.patchPhase;
        });

        baseShell = shellStdenv: StarPUVersion: extraArgs: pkgs.mkShell.override { stdenv = shellStdenv; } (rec {
            buildInputs = with pkgs; [
                pkg-config
                bear # for emacs eglot
                hwloc
                StarPUVersion
                
                # project libs
                criterion

                # for bash to work properlly inside vscode
                bashInteractive
                moreutils
                # for fuser
                psmisc
                # debug tools
                gdb
                valgrind
                # scripting
                python313
                python313Packages.numpy

                #madagascar.packages.${system}.default
            ] ++ (with cudaPacks; [ 
                cuda_nvcc
                cuda_cudart
                cuda_cccl
                cuda_nvml_dev
                libcublas 
                libcusparse
                libcusolver
                libcufft
            ]);
            # export StarPU and hwloc store locations 
            # for use in vscode intellisence
            GCC_STORE_PATH =         "${pkgs.gcc}/include/";
            STARPU_STORE_PATH =      "${StarPUVersion}/include/starpu/1.4";
            CRITERION_STORE_PATH =   "${pkgs.criterion.dev}/include/";
            HWLOC_STORE_PATH =       "${pkgs.hwloc.dev}/include";
            # cudas
            CUDART_STORE_PATH =      "${pkgs.lib.getInclude cudaPacks.cuda_cudart}/include/";
            NVCC_STORE_PATH =        "${cudaPacks.cuda_nvcc}/include/";
            NVML_STORE_PATH =        "${pkgs.lib.getInclude cudaPacks.cuda_nvml_dev}/include/";
            LIBCUSPARSE_STORE_PATH = "${pkgs.lib.getInclude cudaPacks.libcusparse}/include/";
        } // extraArgs);

        star-fletcher = pkgs.callPackage ./star-fletcher.nix {
            StarPU = starpu-base;
        };

        starpu-cuda = starpu-base.override {
            	enableCUDA = true;
	        compileAsRelease = true;
	        enableTrace = false;
            	maxBuffers = 56;
	        cudaPackages = cudaPacks;
	        stdenv = cudaStdenv;
	    };

        star-fletcher-cuda = pkgs.callPackage ./star-fletcher.nix {
            StarPU = starpu-base; #starpu-cuda;
            cudaPackages = cudaPacks;
            enableCUDA = true;
            enableTrace = false;
            compileAsRelease = true;
            stdenv = cudaStdenv;
        };

        star-fletcher-cuda-trace = pkgs.callPackage ./star-fletcher.nix {
            StarPU = starpu-base; #starpu-cuda;
            cudaPackages = cudaPacks;
            enableCUDA = true;
            enableTrace = true;
            compileAsRelease = true;
            stdenv = cudaStdenv;
        };

        star-fletcher-cuda-no-cpu-kernel = pkgs.callPackage ./star-fletcher.nix {
            StarPU = starpu-base; #starpu-cuda;
            cudaPackages = cudaPacks;
            enableCUDA = true;
            enableTrace = false;
            disableCPUKernel = true;
            compileAsRelease = true;
            stdenv = cudaStdenv;
        };

        nixglhost = nix-gl-host.defaultPackage.${system};

        kernel-test = pkgs.callPackage ./tests/kernel-opt/kernel-test.nix {
	    root = self;
            compileAsRelease = true;
            stdenv = cudaStdenv;
            cudaPackages = cudaPacks;
        };
    in
    {
        devShells = {
          default = baseShell cudaStdenv starpu-cuda {};
          no-cuda = baseShell pkgs.gcc13Stdenv (starpu-base.override {
            	enableCUDA = false;
	            compileAsRelease = true;
	            enableTrace = false;
            	maxBuffers = 56;
	        }) {};
            /*
          eztrace-test = pkgs.mkShell {
                buildInputs = [
                    (eztrace.packages.${system}.new-starpu.override {
                        enableTrace = false;
                        maxBuffers = 56;
                        enableCUDA = false;
                        compileAsRelease = false;
                        stdenv = pkgs.gcc13Stdenv;
                    })
                    (eztrace.packages.${system}.eztrace.override {
                        pallas = eztrace.packages.${system}.pallas.override {
                            stdenv = pkgs.gcc13Stdenv;
                        };
                        stdenv = pkgs.gcc13Stdenv;
                        compileAsDebug = true;
                    })
                    pkgs.pkg-config
                    pkgs.hwloc
                    pkgs.python313
                    pkgs.gdb
                ];
            };*/
            cuda-test = pkgs.mkShell.override { stdenv = cudaStdenv; } {
                buildInputs = [ 
                    pkgs.python313
                    pkgs.pkg-config
                    pkgs.hwloc

                    cudaPacks.cuda_gdb
                    starpu-cuda
                    #pkgs.gdb
                    nixglhost 
                    cudaPacks.cuda_cuobjdump
		    cudaPacks.nsight_compute
                ] ++ (with cudaPacks; [
                    cuda_nvcc
                ]);
            };
            star-compilation = pkgs.mkShell {
                buildInputs = [ 
                    star-fletcher-cuda		
                    nixglhost 
                ];
            };

            star-fletcher-trace = pkgs.mkShell {
                buildInputs = [ 
                    star-fletcher-cuda-trace
                    nixglhost 
                ];
            };
            kernel-test-grace = pkgs.mkShell {
                buildInputs = [ 
                    nixglhost 
		    (kernel-test.overrideAttrs { doCheck = false; })
                ];
            };
        };
        packages = {
          default = star-fletcher;
          inherit star-fletcher star-fletcher-cuda star-fletcher-cuda-no-cpu-kernel star-fletcher-cuda-trace kernel-test;
        };
    });
}
