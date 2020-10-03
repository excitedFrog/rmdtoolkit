# Python 3.6.1

import sys


class ParserBase(object):
    def __init__(self, file,
                 tell_time=False, tell_freq=100,
                 err=sys.stderr, out=sys.stdout):
        """
        :param file: the File object to be parsed.
        """
        self.file = file

        self.line = str()  # current line
        self.frame = list()  # current frame, a list of lines
        self.last_pos = int()  # last offset of file
        self.current_pos = int()  # current offset of file
        self.frame_time = int()

        self.tell_time = tell_time
        self.tell_freq = tell_freq
        self.err = err
        self.out = out
        self.tell_format = '{ClassName} Parsing Frame {Time}'.format(ClassName=self.__class__.__name__,
                                                                     Time='{Time}')

    def read_line(self):
        status = 0
        self.line = self.file.readline()
        self.last_pos = self.current_pos
        self.current_pos = self.file.tell()
        if not self.line:
            status = 1
        return status

    def read_lines(self):
        self.frame = self.file.readlines()

    def read_frame(self, **kwargs):
        raise NotImplementedError("Override this function in each subclass.")

    def parse_frame(self, **kwargs):
        raise NotImplementedError("Override this function in each subclass.")

    def parse_file(self, **kwargs):
        raise NotImplementedError("Override this function in each subclass.")

    def time_tell(self):
        if self.tell_time:
            if self.frame_time % self.tell_freq == 0:
                print(self.tell_format.format(Time=self.frame_time), file=self.out)

    @staticmethod
    def split_list(a_list):
        a_list = list(map(lambda x: x.strip('\n'), a_list))
        a_list = list(map(lambda x: x.split(), a_list))
        return a_list
