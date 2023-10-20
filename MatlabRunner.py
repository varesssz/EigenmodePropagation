from subprocess import run
import os.path


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
            "-r", "run('cd(['%s'])');" % self.current_directory,
        ])

    def run_matlab_script(
            self,
            script_fname: str,
    ):
        run([
            self.program_path, "-nodisplay", "-nosplash", "-nodesktop",
            "-r", "run('%s'); exit;" % os.path.join(self.current_directory, script_fname),
        ])

    def export_fixer(self):
        data_folder = os.path.join(self.current_directory, "data")
        for file_name in os.listdir(data_folder):
            if file_name.find(".txt") != -1:
                # Read the contents of file in a variable
                with open(os.path.join(data_folder, file_name), "r") as file:
                    content = file.read()

                # Replace characters
                content = content.replace("i", "j")
                content = content.replace("+-", "-")

                # Open the file in write mode and write the new fixed content
                with open(os.path.join(data_folder, file_name), "w") as file:
                    file.write(content)
