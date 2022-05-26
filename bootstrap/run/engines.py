import os
import time
import datetime
import dataclasses

import rich

import common

import run.queues   as queues
import run.mpi_bins as mpi_bins

@dataclasses.dataclass
class Engine:
    name: str
    slug: str

    def get_targets(self, targets: list) -> list:
        raise common.MFCException(f"MFCEngine::get_targets: not implemented for {self.name}.")

    def run(self, mfc, target_name: str, mpibin: mpi_bins.MPIBinary) -> None:
        raise common.MFCException(f"MFCEngine::run: not implemented for {self.name}.")

    def validate_job_options(self, mfc) -> None:
        raise common.MFCException(f"MFCEngine::validate_job_options: not implemented for {self.name}.")


class SerialEngine(Engine):
    def __init__(self) -> None:
        super().__init__("Serial", "serial")

    def get_targets(self, targets: list) -> list:
        if targets[0] == "mfc":
            return ["pre_process", "simulation", "post_process"]

        return targets

    def run(self, mfc, target_name: str, mpibin: mpi_bins.MPIBinary) -> None:
        self.mfc = mfc

        date = f"> > [bold cyan][{common.get_datetime_str()}][/bold cyan]"
        rich.print(f"{date} Running...")

        start_time = time.monotonic()
        common.execute_shell_command(self.mfc.run.get_exec_cmd(target_name, mpibin))
        end_time   = time.monotonic()

        rich.print(f"> > Done [bold green]✓[/bold green] (in {datetime.timedelta(seconds=end_time - start_time)})")

    def validate_job_options(self, mfc) -> None:
        if mfc.args["nodes"] != 1:
            raise common.MFCException("SerialEngine: In serial mode, only node can be used.")


class ParallelEngine(Engine):
    def __init__(self) -> None:
        super().__init__("Parallel", "parallel")

    def get_targets(self, targets: list) -> list:
        if targets[0] == "mfc" or len(targets) != 1:
            raise common.MFCException("The parallel engine requires a unique target to run.")

        return targets

    def run(self, mfc, target_name: str, mpibin: mpi_bins.MPIBinary) -> None:
        self.mfc = mfc

        system = queues.get_system()
        rich.print(f"> > Detected the [bold magenta]{system.name}[/bold magenta] queue system.")

        self.create_batch_file(system, target_name, mpibin)

        self.execute_batch_file(system, target_name)

    def get_batch_filepath(self, system: queues.QueueSystem, target_name: str):
        case_dirpath = self.mfc.run.get_case_dirpath()

        return os.path.abspath(f"{case_dirpath}/{target_name}.sh")

    def create_batch_file(self, system: queues.QueueSystem, target_name: str, mpibin: mpi_bins.MPIBinary):
        job_name = self.mfc.run.get_job_name(target_name)

        load_modules: str = "module load mpi"

#        if common.does_cmd_exist("module"):
#            modules: list = common.loaded_modules()
#            load_modules: str = f"""\
#echo " :) Loading modules... ({common.format_list_to_string(modules)})"
#
#module purge || true
#{chr(10).join([ f"module load {m} || true" for m in modules ])}
#
#echo ""
#"""

        BATCH_CONTENT: str = f"""\
#!/usr/bin/env bash
{system.gen_batch_header(self.mfc.args, job_name)}

TABLE_FORMAT_LINE="| - %-14s %-35s - %-14s %-35s |\\n"
TABLE_HEADER="+-----------------------------------------------------------------------------------------------------------+ \\n"
TABLE_FOOTER="+-----------------------------------------------------------------------------------------------------------+ \\n"
TABLE_TITLE_FORMAT="| %8s %-96s |\\n"
TABLE_CONTENT=$(cat <<-END
$(printf "$TABLE_FORMAT_LINE" "Start-time:"    "$(date +%T)"                       "Start-date:"    "$(date +%T)")
$(printf "$TABLE_FORMAT_LINE" "Partition:"     "{self.mfc.args["partition"]}"      "Walltime:"      "{self.mfc.args["walltime"]}")
$(printf "$TABLE_FORMAT_LINE" "Account:"       "{self.mfc.args["account"]}"        "Nodes:"         "{self.mfc.args["nodes"]}")
$(printf "$TABLE_FORMAT_LINE" "CPUs (/node):"  "{self.mfc.args["cpus_per_node"]}"  "GPUs (/node):"  "{self.mfc.args["gpus_per_node"]}")
$(printf "$TABLE_FORMAT_LINE" "Input File:"    "{self.mfc.args["input"]}"          "Engine"         "{self.mfc.args["engine"]}")
$(printf "$TABLE_FORMAT_LINE" "Queue System:"  "{system.name}"                     "Mode:"          "{self.mfc.args["mode"]}")
$(printf "$TABLE_FORMAT_LINE" "Email:"         "{self.mfc.args["email"]}"          "Job Name:"      "{job_name}")
$(printf "$TABLE_FORMAT_LINE" "MPI Binary:"    "{mpibin.bin}"                      ""      "")\\n
END
)

printf "$TABLE_HEADER"
printf "$TABLE_TITLE_FORMAT" "Starting" "{job_name}:"
printf "$TABLE_CONTENT"
printf "$TABLE_FOOTER\\n"

t_start=$(date +%s)

{load_modules}\

echo -e ":) Running MFC..."
{self.mfc.run.get_exec_cmd(target_name, mpibin)}

code=$?

t_stop="$(date +%s)"

status_msg="SUCCESS"
if [ "$code" -ne "0" ]; then
    status_msg="FAILURE"
fi

printf "\\n$TABLE_HEADER"
printf "$TABLE_TITLE_FORMAT" "Finished" "{job_name}:"
printf "$TABLE_FORMAT_LINE" "Total-time:"  "$(expr $t_stop - $t_start)s"  "Status:"    "$status_msg"
printf "$TABLE_FORMAT_LINE" "End-time:"    "$(date +%T)"                  "End-date:"  "$(date +%T)"
printf "$TABLE_FOOTER"

printf "\\nThank you for using MFC !\\n\\n"

exit $code
"""

        filepath = self.get_batch_filepath(system, target_name)

        common.file_write(filepath, BATCH_CONTENT)

    def execute_batch_file(self, system: queues.QueueSystem, target_name: str):
        if 0 != os.system(system.gen_submit_cmd(self.get_batch_filepath(system, target_name))):
            raise common.MFCException(f"Running batch file for {system.name} failed.")

    def validate_job_options(self, mfc) -> None:
        pass


ENGINES = [ SerialEngine(), ParallelEngine() ]

def get_engine(slug: str) -> Engine:
    engine: Engine = None
    for candidate in ENGINES:
        candidate: Engine

        if candidate.slug == slug:
            engine = candidate
            break

    if engine == None:
        raise common.MFCException(f"Unsupported engine {slug}.")

    return engine