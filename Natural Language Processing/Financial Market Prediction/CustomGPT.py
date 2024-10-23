class CustomGPT(nn.Module):
  '''
  This Model has same architecture as GPT-2 from Huggingface.
  However here we use two (encoder, decoder) layers to create token embeddings (768 dimension).

  '''
  def __init__(self, gpt_model, window_size):
    super().__init__()
    self.decoder = nn.Sequential(
        nn.Linear(4, 256),
        nn.BatchNorm1d(window_size, eps=1e-08),
        nn.ReLU(),
        nn.Linear(256, 512),
        nn.BatchNorm1d(window_size, eps=1e-08),
        nn.ReLU(),
        nn.Linear(512, 768)
    )
    self.gpt = gpt_model
    self.encoder = nn.Sequential(
        nn.Linear(768, 256),
        nn.BatchNorm1d(window_size, eps=1e-08),
        nn.ReLU(),
        nn.Linear(256, 128),
        nn.BatchNorm1d(window_size, eps=1e-08),
        nn.ReLU(),
        nn.Linear(128, 4)
    )

  def forward(self, inputs_embeds):
    x = self.decoder(inputs_embeds)
    x = self.gpt(inputs_embeds=x).last_hidden_state
    x = self.encoder(x)
    return x