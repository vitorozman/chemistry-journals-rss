
import re

def eq_name(res_name, auth_name):
    res_names = re.split(r'[-. ]', res_name.lower())
    auth_names = re.split(r'[-. ]', auth_name.lower())
    res_names = [n for n in res_names if n.strip()]
    auth_names = [n for n in auth_names if n.strip()]
    if len(auth_names)==0 or len(res_names)==0:
        return False
    for res in res_names.copy():
        for j, auth in enumerate(auth_names):
            k = min(len(res), len(auth))
            if k > 1:
                if res == auth:
                    auth_names.pop(j)
                    res_names.pop(res_names.index(res))
                    break
            elif res[:k] == auth[:k]:
                auth_names.pop(j)
                res_names.pop(res_names.index(res))
                break
    if len(res_names) == 0 or len(auth_names) == 0:
        return True
    return False



def normalize_institution(text: str) -> str:
    stop_words = {
        "the", "of", "and", "for", "in", "at", "on", "to", "a", "an"
    }
    if not text:
        return ""
    text = text.lower()
    text = re.sub(r"[^\w\s]", " ", text)
    tokens = [t for t in text.split() if t not in stop_words]
    return " ".join(tokens)

def eq_affiliation(res_inst:str, auth_inst:str) -> bool:
    """This function is used when comparing researcher institution in db and affiliation from outsource"""
    res_insts = re.split(r'[- ,.]', res_inst.lower())
    auth_insts = re.split(r'[- ,.]', auth_inst.lower())
    for inst_part in res_insts.copy():
        if inst_part not in auth_insts:
            return False
    return True


def eq_institution(db_inst:str, new_inst:str) -> bool:
    """This function is used when comparing institution in db and new institution"""
    norm_db_inst = normalize_institution(db_inst)
    norm_new_inst = normalize_institution(new_inst)
    return norm_db_inst == norm_new_inst


