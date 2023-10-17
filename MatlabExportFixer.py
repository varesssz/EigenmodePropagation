import os.path

if __name__ == '__main__':
    folder = "data/"
    for file_name in os.listdir(folder):
        if file_name.find(".txt") != -1:
            # Read the contents of file in a variable
            with open(os.path.join(folder, file_name), "r") as file:
                content = file.read()

            # Replace characters
            content = content.replace("i", "j")

            # Open the file in write mode and write the new fixed content
            with open(os.path.join(folder, file_name), "w") as file:
                file.write(content)
