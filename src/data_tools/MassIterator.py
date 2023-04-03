
class MassIterator:
    def __init__(self, df):
        self.df = df
        self.masses = df[df.type=='sig'].mass.unique()
        self.current = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self.current < len(self.masses):
            mass = self.masses[self.current]
            self.current += 1
            return mass, self.df[self.df.mass==mass]
        self.current=0
        raise StopIteration
