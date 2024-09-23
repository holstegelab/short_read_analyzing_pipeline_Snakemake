import sys
import gzip
import os
import threading
import time
      
class PairedFastQReader(threading.Thread):
    def __init__(self, filea,fileb):
        self.filea = filea
        self.fileb = fileb
        self.buffer = {}
        self.splitchar = None
        self.buffer_size = 1000

        self.file_finished = False
        self.continue_reading = threading.Event()
        self.search_read = None
        self.read_found_or_buffer_full = threading.Event()

        threading.Thread.__init__(self)


    def run(self):
        with os.popen('pigz -dc ' + self.filea) as f1:
            with os.popen('pigz -dc ' + self.fileb) as f2:
                try:
                    f1 = iter(f1)
                    f2 = iter(f2)
                    while True:
                        header1 = f1.__next__().strip()
                        seq1 = f1.__next__().strip()
                        dummy = f1.__next__()
                        qual1 = f1.__next__().strip()
                        header2 = f2.__next__().strip()
                        seq2 = f2.__next__().strip()
                        dummy = f2.__next__()
                        qual2 = f2.__next__().strip()
                        assert header1.startswith('@'), 'Parsing error'
                        assert header2.startswith('@'), 'Parsing error'
                        if self.splitchar is None:
                            if ' ' in header1:
                                self.splitchar = ' '
                            elif '/' in header1:
                                self.splitchar = '/'
                            else:
                                raise RuntimeError('Cannot determine FastQ header split character')
                                
                        name1, r1 = header1.split(self.splitchar)
                        name2, r2 = header2.split(self.splitchar)
                        name1 = name1[1:] #strip '@'
                        name2 = name2[1:] #strip '@'

                        assert name1 == name2, 'FastQ input files not in same read order'

                        self.buffer[name1] = (name1, seq1, seq2, qual1, qual2)
                        
                        if name1 == self.search_read or self.search_read is True:
                            self.read_found_or_buffer_full.set()

                        if len(self.buffer) >= self.buffer_size:
                            self.read_found_or_buffer_full.set()
                            self.continue_reading.clear()
                            self.continue_reading.wait()

                except StopIteration:
                    self.read_found_or_buffer_full.set()
                    self.file_finished = True

    def popRead(self, name):
        while not name in self.buffer:
            self.search_read = name
            self.read_found_or_buffer_full.clear()
            while True:
                self.read_found_or_buffer_full.wait(timeout=5)
                if name in self.buffer:
                    break
                elif self.file_finished:
                    raise RuntimeError(f"Read {name} not found while file has ended. Is BAM file reordered?")
                elif len(self.buffer) >= self.buffer_size: #buffer full, read not found
                    self.buffer_size *=2
                    sys.stderr.write(f"Order of reads has changed between BAM and FASTQ. Growing look-back buffer to {self.buffer_size}\n")
                    sys.stderr.flush()
                    self.continue_reading.set()
                    break
                else:    
                    self.continue_reading.set()
                    #wait a bit more
                    continue

        read = self.buffer[name]
        del self.buffer[name]
        
        if len(self.buffer) < self.buffer_size:
            self.continue_reading.set()
        return read

    def retrieveRead(self):
        if len(self.buffer) < self.buffer_size:
            self.continue_reading.set()

        if not self.buffer:
            if self.file_finished:
                raise StopIteration()
                
            self.read_found_or_buffer_full.clear()
            self.search_read = True
            self.read_found_or_buffer_full.wait()
            if not self.buffer:
                assert self.file_finished
                raise StopIteration()

        return self.buffer.popitem()[1]


class PairedFastQReaderSimple:
    def __init__(self, filea,fileb):
        self.filea = filea
        self.fileb = fileb
        self.splitchar = None

        self.file_finished = False
        self.search_read = None


    def retrieveRead(self):
        with os.popen('pigz -dc ' + self.filea) as f1:
            with os.popen('pigz -dc ' + self.fileb) as f2:
                try:
                    f1 = iter(f1)
                    f2 = iter(f2)
                    while True:
                        header1 = f1.__next__().strip()
                        seq1 = f1.__next__().strip()
                        dummy = f1.__next__()
                        qual1 = f1.__next__().strip()
                        header2 = f2.__next__().strip()
                        seq2 = f2.__next__().strip()
                        dummy = f2.__next__()
                        qual2 = f2.__next__().strip()
                        assert header1.startswith('@'), 'Parsing error'
                        assert header2.startswith('@'), 'Parsing error'
                        if self.splitchar is None:
                            if ' ' in header1:
                                self.splitchar = ' '
                            elif '/' in header1:
                                self.splitchar = '/'
                            else:
                                raise RuntimeError('Cannot determine FastQ header split character')
                                
                        name1, r1 = header1.split(self.splitchar)
                        name2, r2 = header2.split(self.splitchar)
                        name1 = name1[1:] #strip '@'
                        name2 = name2[1:] #strip '@'

                        assert name1 == name2, 'FastQ input files not in same read order'

                        yield (name1, seq1, seq2, qual1, qual2)

                except StopIteration:
                    pass


