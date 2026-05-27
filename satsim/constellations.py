from pkg_resources import resource_filename

# Define the path to the data directory
data_path = resource_filename("satsim", "data/")

def filter_dtc_satellites(input_file, output_file):
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            lines = infile.readlines()
            for i in range(0, len(lines), 3):
                if "[DTC]" in lines[i]:
                    outfile.write(lines[i])
                    outfile.write(lines[i + 1])
                    outfile.write(lines[i + 2])

    except FileNotFoundError:
        print(f"Error: The file '{input_file}' was not found.")
    except IndexError:
        print("Error: The file does not have complete TLE data for all satellites.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

def filter_oneweb_satellites(input_file, output_file, satellites_excluded):
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            lines = infile.readlines()
            for i in range(0, len(lines), 3):
                satellite_number = lines[i].split()[0]
                if satellite_number not in satellites_excluded:
                    outfile.write(lines[i])
                    outfile.write(lines[i + 1])
                    outfile.write(lines[i + 2])

    except FileNotFoundError:
        print(f"Error: The file '{input_file}' was not found.")
    except IndexError:
        print("Error: The file does not have complete TLE data for all satellites.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

input_file = "{}starlink.txt".format(data_path)
output_file = "{}starlink_dtc.txt".format(data_path)
#filter_dtc_satellites(input_file, output_file)

input_file = "{}oneweb_tle.txt".format(data_path)
output_file = "{}oneweb.txt".format(data_path)
satellites_excluded = ["ONEWEB-0050", "ONEWEB-0693", "ONEWEB-0687", "ONEWEB-0688", "ONEWEB-0015", "ONEWEB-0690", "ONEWEB-0691", "ONEWEB-0695", "ONEWEB-0689", "ONEWEB-0702", "ONEWEB-0692", "ONEWEB-0698", "ONEWEB-0697", "ONEWEB-0618", "ONEWEB-0013", "ONEWEB-0699", "ONEWEB-0700", "ONEWEB-0701", "ONEWEB-0704", "ONEWEB-0705", "ONEWEB-0706", "ONEWEB-0707", "ONEWEB-0708", "ONEWEB-0721"]
#filter_oneweb_satellites(input_file, output_file, satellites_excluded)
