dry_run = False
project = "K6-ALL"
checkpoint_dir = f"checkpoints/{project}"

hours_per_eval = 4
seconds_per_step = 2
steps_per_eval = (hours_per_eval * 60 * 60) * seconds_per_step

model = BertForSequenceClassification.from_pretrained(model_path, num_labels=4)
training_args = TrainingArguments(
    output_dir=checkpoint_dir,
    num_train_epochs=10,
    #evaluation_strategy="epoch",
    evaluation_strategy="steps",
    eval_steps=steps_per_eval,
    #max_steps=100,
    report_to="wandb",
    per_device_train_batch_size=15,
    per_device_eval_batch_size=128,
    logging_strategy="steps",
    logging_first_step=True,
    logging_steps=1,
)

# https://scikit-learn.org/stable/modules/generated/sklearn.metrics.f1_score.html

metric_acc = load_metric("accuracy")
metric_f1 = load_metric("f1")

def compute_metrics(eval_pred):
    logits, labels = eval_pred
    predictions = np.argmax(logits, axis=-1)
    return dict(
        Accuracy=metric_acc.compute(predictions=predictions, references=labels),
        F1=metric_f1.compute(predictions=predictions, references=labels, average="macro"),
    )

trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=ds_train,
    eval_dataset=ds_eval,
    compute_metrics=compute_metrics,
)

if dry_run:
    os.environ["WANDB_MODE"] = "offline"
    os.environ["WANDB_DISABLED"] = "true"
else:
    if "WANDB_MODE" in os.environ:
        del os.environ["WANDB_MODE"]
    if "WANDB_DISABLED" in os.environ:
        del os.environ["WANDB_DISABLED"]

try:
    with wandb.init() as run:
        res = trainer.train()
        print(res)
finally:
    if not dry_run:
        path_final = f"{checkpoint_dir}/final"
        model.save_pretrained(path_final)
    del model
