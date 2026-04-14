import Pluto
import PlutoSliderServer

# 1. Configuration: Control the processing load
pluto_config = Pluto.Configuration.from_flat_kwargs(
    capability_maximum_active_notebooks = 2, # ONLY process 2 notebooks at once
    workspace_use_distributed = true
)

# 2. Server settings
server_options = PlutoSliderServer.Configuration.ServerOptions(
    port = 8081,
    host = "0.0.0.0"
)

# 3. Launch from the Git URL
# This will clone the repo into a temporary subfolder in /home/user/pluto-server/
PlutoSliderServer.run_git_directory(
    "https://github.com/your-username/your-notebook-repo.git";
    pluto_configuration = pluto_config,
    slider_server_configuration = server_options,
    static_export = false
)
