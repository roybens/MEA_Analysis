import tkinter as tk
from tkinter import filedialog
import json

# This class creates a GUI for selecting multiple directories
class MultiDirPicker(tk.Frame):
    def __init__(self, master=None):
        super().__init__(master)
        self.master = master
        self.pack()
        self.create_widgets()

    # This method creates all the widgets for the GUI
    def create_widgets(self):
        # Button for selecting directories
        self.select_button = tk.Button(self, font=("Helvetica", 16))
        self.select_button["text"] = "Select Directories"
        self.select_button["command"] = self.select_directories
        self.select_button.pack(side="top")

        # Listbox for displaying selected directories
        self.directory_list = tk.Listbox(self, width=100, font=("Helvetica", 16))
        self.directory_list.pack(side="top")

        # Button for removing selected directory
        self.remove_button = tk.Button(self, font=("Helvetica", 16))
        self.remove_button["text"] = "Remove Selected Directory"
        self.remove_button["command"] = self.remove_directory
        self.remove_button.pack(side="top")

        # Button for closing the GUI
        self.quit = tk.Button(self, text="OK", fg="red", font=("Helvetica", 16),
                              command=self.master.destroy)
        self.quit.pack(side="bottom")

    # This method opens a directory picker and adds the selected directory to the list
    def select_directories(self):
        initial_dir = "/mnt" if not self.master.selected_directories else None
        folder_path = filedialog.askdirectory(initialdir=initial_dir)
        if folder_path:
            self.master.selected_directories.append(folder_path)
            self.directory_list.insert(tk.END, folder_path)

    # This method removes the selected directory from the list
    def remove_directory(self):
        selected_indices = self.directory_list.curselection()
        for index in reversed(selected_indices):
            self.master.selected_directories.pop(index)
            self.directory_list.delete(index)

# This function creates the GUI and returns the selected directories
def select_folders():
    root = tk.Tk()
    root.title("Roy Ben-Shalom Lab Axon Reconstruction Pipeline")
    root.selected_directories = []
    app = MultiDirPicker(master=root)

    # Center the window
    window_width = 800
    window_height = 400
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()
    position_top = int(screen_height / 2 - window_height / 2)
    position_right = int(screen_width / 2 - window_width / 2)
    root.geometry(f"{window_width}x{window_height}+{position_right}+{position_top}")

    app.mainloop()
    return root.selected_directories

# This function prints the selected directories and returns them as a JSON string
def main(pre_selected_folders = None, debug_mode = False):    
    if debug_mode:
        debug_folders = pre_selected_folders
        print("\n\033[1;33m*** DEBUG MODE ENABLED in file_selector.py ***\033[0m") 
        for folder in debug_folders:
            print(folder)
        selected_folders = debug_folders
    else:
        selected_folders = select_folders()        
        for folder in selected_folders:
            print(folder)
    print(f"\nFolders selected for analysis:")
    return selected_folders

# If this script is run directly, call the main function
if __name__ == "__main__":
    main()