from slack_sdk import WebClient
from slack_sdk.errors import SlackApiError

def send_message_to_slack(text, channel):
    client = WebClient(token='xoxb-1294950085457-2100256556481-ylSxWRn0YPV1itPoSa99XJHg')

    try:
        response = client.chat_postMessage(channel='#' + channel, text=text)
        assert response["message"]["text"] == text
    except SlackApiError as e:
        # You will get a SlackApiError if "ok" is False
        assert e.response["ok"] is False
        assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'
        print(f"Got an error: {e.response['error']}")

def send_file_to_slack(filepath, channel):

    client = WebClient(token='xoxb-1294950085457-2100256556481-ylSxWRn0YPV1itPoSa99XJHg')

    try:
        #filepath="./tmp.txt"
        response = client.files_upload(channels='#' + channel, file=filepath)
        assert response["file"]  # the uploaded file
    except SlackApiError as e:
        # You will get a SlackApiError if "ok" is False
        assert e.response["ok"] is False
        assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'
        print(f"Got an error: {e.response['error']}")