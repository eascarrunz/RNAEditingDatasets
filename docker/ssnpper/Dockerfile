FROM julia:1.10.9-alpine3.21
WORKDIR /work

COPY ssnpper.jl /script/ssnpper.jl
COPY Project.toml /script/Project.toml
COPY docker/ssnpper/project_environment.jl /script/project_environment.jl
RUN cd /script && julia -q --project="@." project_environment.jl
COPY --chmod=555 docker/ssnpper/launch.sh /script/launch.sh
RUN ln -n /script/launch.sh /usr/bin/ssnpper

CMD ["sh"]
