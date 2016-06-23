# standard libraries
import io

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"


class BasicFastq:
    def __init__(self):
        self.clear()

    def add_line(self, new_line):
        self.lines.append(new_line)
        result = True if len(self.lines) >= 4 else False
        return result

    @property
    def sequence(self):
        return self.lines[1].strip()

    @sequence.setter
    def sequence(self, value):
        self.lines[1] = value + "\n"


    @property
    def quality(self):
        return self.lines[3].strip()

    @quality.setter
    def quality(self, value):
        self.lines[3] = value + "\n"

    def to_string(self):
        return "".join(self.lines)

    def clear(self):
        self.lines = []


class FastqHandler:
    file_suffix = "fastq"

    def __init__(self, file_path, input_is_fastq_string=False):
        if not input_is_fastq_string:
            self.filepath = file_path
            self.generator = open(self.filepath, 'r')
        else:
            self.filepath = None
            self.generator = io.StringIO(file_path)

        self.basic_fastq = BasicFastq()
        self.is_done = False

    def get_next(self):
        curr_line = next(self.generator)
        self.is_done = self.basic_fastq.add_line(curr_line)

    def clear_record(self):
        self.basic_fastq.clear()

    def close(self):
        self.generator.close()


def paired_fastq_generator(fw_fastq_handler, rv_fastq_handler, get_full_record=False):
    fastq_handlers = [fw_fastq_handler, rv_fastq_handler]

    while True:
        try:
            for curr_handler in fastq_handlers:
                curr_handler.get_next()
        except StopIteration:
            for curr_handler in fastq_handlers:
                if not curr_handler.is_done:
                    print(("Error: {0} generator terminated"
                           "with unfinished record").format(
                        curr_handler.filepath))
            break

        done_booleans = [x.is_done for x in fastq_handlers]
        if True in done_booleans:
            # at least one handler thinks its record is done
            if False in done_booleans:
                # at least one handler DOESN'T think its record is done
                print("Error--fastq records aren't in sync")
                break
            else:
                if get_full_record:
                    result = (fw_fastq_handler.basic_fastq, rv_fastq_handler.basic_fastq)
                else:
                    result = (fw_fastq_handler.basic_fastq.sequence.upper(),
                              rv_fastq_handler.basic_fastq.sequence.upper())
                yield result

                # prepare for next fastq record
                for curr_handler in fastq_handlers:
                    curr_handler.clear_record()

    for curr_handler in fastq_handlers:
        curr_handler.close()