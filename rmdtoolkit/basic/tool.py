# Python 3.6.1

import os
import time
import subprocess
from itertools import zip_longest

import numpy as np
from twilio.rest import Client

from rmdtoolkit.config.config_manager import ConfigManager

CM = ConfigManager()

TWILIO_SID = CM.get('twilio', 'twilio_sid')
TWILIO_AUTH_TOKEN = CM.get('twilio', 'twilio_token')
TWILIO_NUMBER = CM.get('twilio', 'twilio_virtual_number')

QOS = CM.get('slurm', 'special_partition_prefix')
PARTITION_LIST = CM.get_list('slurm', 'partition_list')
DEFAULT_PARTITION = CM.get('slurm', 'default_partition')


def find_file(extension_name, directory='./'):
    file_list = os.listdir(directory)
    return [_.split('.')[0] for _ in file_list if _.endswith(extension_name)]


def list_filter(full_list, exclude_list):
    return [x for x in full_list if x not in exclude_list]


def string_is_true(s):
    if s.lower() in ['yes', 'true']:
        return True
    else:
        return False


def second_largest(numbers):
    count = 0
    m1 = m2 = float('-inf')
    for x in numbers:
        count += 1
        if x > m2:
            if x >= m1:
                m1, m2 = x, m1
            else:
                m2 = x
    return m2 if count >= 2 else None


# Get list average for lists with different lengths.
def list_average(lst_of_lst):
    average = [np.ma.average(np.ma.masked_values(temp_list, None)) for temp_list in zip_longest(*lst_of_lst)]
    return average


def raise_time():
    time_info = time.localtime()
    return '%s/%s/%s %s:%s:%s' \
           % (time_info[1], time_info[2], time_info[0],
              str(time_info[3]).rjust(2, '0'), str(time_info[4]).rjust(2, '0'), str(time_info[5]).rjust(2, '0'))


def log_error(error_origin, error_msg, error_lv=''):
    with open('%szerror.log' % CM.get('basics', 'error_save_dir'), 'a') as error_file:
        error_file.write('**********\n[%s ERROR] @ %s\nFrom %s:\n %s'
                         % (error_lv, raise_time(), error_origin, error_msg))


def pbc_dist(pos1, pos2, boxlen, retarray=False):
    r = pos1 - pos2
    r = np.array([abs(ri) if abs(ri) < boxlen[i]/2 else boxlen[i] - abs(ri) for i, ri in enumerate(r)])
    if retarray:
        return r
    else:
        dist = np.sqrt(np.dot(r, r))
        return dist


def check_pid(pid):
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    else:
        return True


def text_inform(phone_num, message):
    client = Client(TWILIO_SID, TWILIO_AUTH_TOKEN)
    client.api.account.messages.create(to=phone_num, from_=TWILIO_NUMBER, body=message)


def law_of_cosines(oa, ob):
    ab = np.squeeze(np.asarray(oa - ob))
    ab2 = np.dot(ab, ab)
    oa = np.squeeze(np.asarray(oa))
    ob = np.squeeze(np.asarray(ob))
    oa2 = np.dot(oa, oa)
    oa = np.sqrt(oa2)
    ob2 = np.dot(ob, ob)
    ob = np.sqrt(ob2)
    theta = np.arccos((oa2 + ob2 - ab2) / (2 * oa * ob))
    return theta


def find_idle_partition():
    idle_partition = str()
    idle_amount = 0
    for partition in PARTITION_LIST:
        try:
            info = subprocess.check_output('sinfo -p %s|grep idle' % partition, shell=True)
        except subprocess.CalledProcessError:
            continue
        info = info.split()
        if int(info[3].decode()) > idle_amount:
            idle_amount = int(info[3].decode())
            idle_partition = info[0].decode()
    if idle_partition == str():
        idle_partition = DEFAULT_PARTITION
    idle_partition = idle_partition.strip('*')
    return idle_partition


def chunks(l, n):
    n = max(1, n)
    return [l[i:i+n] for i in range(0, len(l), n)]
