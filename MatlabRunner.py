from subprocess import run
import os.path
import time


class MatlabRunner:
    def __init__(
            self,
            program_path="/home/varesz/MATLAB/bin/matlab",
    ):
        self.program_path = program_path
        self.current_directory = os.path.join(os.getcwd(), "MATLAB")

    def open_matlab(self):
        run([
            self.program_path, "-nosplash", "-desktop",
            "-r", "cd(['%s'])" % self.current_directory,
        ])

    def run_matlab_script(
            self,
            script_fname: str,
    ):
        start_sec = time.time()
        run([
            self.program_path, "-nodisplay", "-nosplash", "-nodesktop",
            "-r", "run('%s'); exit;" % os.path.join(self.current_directory, script_fname),
        ])
        elapsed_sec = time.time() - start_sec
        print("=======================================")
        print(" The script took %dh %dmin %dsec to run" % (elapsed_sec // 3600, elapsed_sec // 60, elapsed_sec % 60))
        print("=======================================\n\n")

    def export_folder_fixer(
            self,
            sub_folder_in_cd="data",
    ):
        for file_name in os.listdir(os.path.join(self.current_directory, sub_folder_in_cd)):
            if file_name.find(".txt") != -1:
                self.export_fixer(
                    file_path_in_cd=os.path.join(sub_folder_in_cd, file_name)
                )

    def export_fixer(
            self,
            file_path_in_cd: str,
    ):
        file_path = os.path.join(self.current_directory, file_path_in_cd)
        # Read the contents of file in a variable
        with open(file_path, "r") as file:
            content = file.read()

        # Replace characters
        content = content.replace("i", "j")
        content = content.replace("+-", "-")

        # Open the file in write mode and write the new fixed content
        with open(file_path, "w") as file:
            file.write(content)
