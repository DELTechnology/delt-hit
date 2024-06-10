from pathlib import Path

from test_counts import read_txt, read_json


def main():

    path = Path('../../data/simulation')

    counts_true = read_txt(path / 'counts_true.txt')[1:]
    counts_pred = read_txt(path / 'counts.txt')[1:]

    structure = read_json(path / 'structure.json')
    codes1 = structure['B1']['Sequences']
    codes2 = structure['B2']['Sequences']

    reads = read_txt(path / 'simulation.fastq')
    # print(len(reads))
    # exit()

    errors = []
    num_reads = len(counts_pred)

    for i, (count_true, count_pred) in enumerate(zip(counts_true, counts_pred)):
        if count_true != count_pred:
            code1_true, code2_true = count_true.split('\t')[1:]
            code1_pred, code2_pred = count_pred.split('\t')[1:]
            errors += [i]
            print(count_true.strip())
            print(count_pred.strip())
            print('True:', codes1[int(code1_true)], codes1[int(code2_true)])
            print('Pred:', codes1[int(code1_pred)], codes1[int(code2_pred)])
            print('\n')


    print(f'{num_reads - len(errors)} / {num_reads}')
    print(errors)


if __name__ == '__main__':

    main()