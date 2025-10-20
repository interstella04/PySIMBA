from .result import Result
from .tool import Tools


def main():
    
    print("Welcome to pysimba! \n You have the following options:")

    print("'1' - Run the fit")
    print("'2' - Add a new measurement")

    choice = input("Enter your choice: ").strip()


    if choice == "1":
        # Runs the Fit. This file gets executed with calling python3 pysimba
        Result.Run()
    elif choice == "2":
        print("Note: To add a new measurement you will need access to the directory 'data/add/'. There you will have to put your pickle file containing the dictionary you want to add. If that is not the case, please do so.")
        added = input("Is your file in 'data/add/' ? (y/n)").strip()
        if added == "y":
            Tools.AddNewMeasurement()
        elif added == "n":
            print("You can't add a new dictionary.")
            return
        else:
            print("Not a valid option. Please use 'y' or 'n'")
    return
    


if __name__ == "__main__":
    main()