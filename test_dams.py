import urllib.request
import json
try:
    with urllib.request.urlopen("http://127.0.0.1:8000/api/dams") as response:
        data = json.loads(response.read().decode())
        for dam in data:
            if dam["alert"] != "NORMAL":
                print(f"ALERT! {dam['name']}: {dam['alert']}, Level: {dam['current_level']}")
        print(f"Total dams fetched: {len(data)}")
        if len(data) > 0:
            print(f"Example dam: {data[0]}")
except Exception as e:
    print(f"Error: {e}")
