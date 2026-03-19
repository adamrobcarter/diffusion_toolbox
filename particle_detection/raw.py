"""
Print the particles array into the console
"""
import common

width = 10
num = 100

if __name__ == '__main__':
    parser = common.argparser()
    args = parser.parse_args()

    for file in args.files:
        data = common.load(f'particle_detection/data/particles_{file}.npz')
        particles = data['particles']
        labels = data['particles_labels']

        for l in labels:
            print(f'{l:>{width}}', end=' ')

        for i in range(num + 1):
            print()
            for j in range(len(labels)):
                print(f'{particles[i, j]:>{width}.4f}', end=' ')
        print()
