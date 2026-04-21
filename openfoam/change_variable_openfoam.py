def change_variable_in_openfoam_file(filename, old_value, new_value):
    with open(filename, "r", encoding="utf-8") as f:
        content = f.read()

    content = content.replace(str(old_value), str(new_value))

    with open(filename, "w") as f:
        f.write(content)
