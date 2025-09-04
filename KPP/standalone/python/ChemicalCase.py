class ChemicalCase:
    def __init__(self, filename):
        self.filename = filename
        self.parse_file()

    def parse_file(self):
        with open(self.filename, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if 'Meteorological and general grid cell metadata' in line:
                    self.parse_meteorological_fields(lines[i+1:i+12])
                elif 'KPP Integrator-specific parameters' in line:
                    self.parse_integrator_specific_parameters(lines[i+1:i+5])
                elif 'CSV data of full chemical state' in line:
                    self.parse_chemical_state(lines[i+6:])

    def parse_meteorological_fields(self, lines):
        for line in lines:
            key, value = line.split(':', 1)
            key = key.strip().replace(' ', '_').replace('(', '').replace(')', '').lower()
            value = self.try_parse_number(value.strip())
            setattr(self, key, value)

    def parse_integrator_specific_parameters(self, lines):
        for line in lines:
            key, value = line.split(':', 1)
            key = key.strip().replace(' ', '_').replace('(', '').replace(')', '').lower()
            value = self.try_parse_number(value.strip())
            setattr(self, key, value)

    def try_parse_number(self, s):
        try:
            return int(s)
        except ValueError:
            try:
                return float(s)
            except ValueError:
                return s

    def parse_chemical_state(self, lines):
        self.species = {}
        self.rate_constants = []
        self.reaction_rates = []
        current_list = None
        for line in lines:
            key, value = line.split(',')
            key = key.strip()
            value = value.strip()
            try:
                value = float(value)
            except ValueError:
                # check if it's a float in scientific notation without the 'e'
                # this can happen for 3 digit exponents
                if '-' in value and '.' in value:
                    mantissa, exponent = value.split('-')
                    try:
                        value = float(mantissa) * 10 ** -int(exponent)
                        print("Warning: scientific notation without 'e' found in file " + self.filename)
                    except ValueError:
                        raise ValueError(f"Error: unable to convert {value} to float in file {self.filename}")
                else:
                    raise ValueError(f"Error: unable to convert {value} to float in file {self.filename}")
            if key == 'R1':
                current_list = self.rate_constants
            elif key == 'A1':
                current_list = self.reaction_rates
            if current_list is None:
                self.species[key] = value
            else:
                current_list.append(value)

# load in a data example: assumed local, though a full path can be specified
# filepath = "Beijing_L1_20200106_1345.txt"
# beijing_surface = ChemicalCase(filepath)